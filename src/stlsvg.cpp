#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Surface_mesh.h>

#include <cmath>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <format>
#include <sstream>
#include <algorithm> // for stable_sort

#include "orient.cpp"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/bind.h>
#endif

#include "clipper/clipper.hpp"

#define LOG(...)                                   \
  do {                                             \
    std::cerr << std::format(__VA_ARGS__) << "\n"; \
  } while (0)

namespace {

using Polyline = std::vector<ClipperLib::DoublePoint>;
using Polylines = std::vector<Polyline>;

struct ViewBox {
  double minX, minY, width, height;
};

bool IsEqual(double a, double b, double tolerance) {
  return std::abs(a - b) <= tolerance;
}

}  // namespace

namespace slice {

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Plane_3 = Kernel::Plane_3;
using Aff_3 = Kernel::Aff_transformation_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using Face_index = Surface_mesh::Face_index;
namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Input>
std::optional<Surface_mesh> ReadSTL(Input& filename_or_buf) {
  Surface_mesh mesh;
  if (!CGAL::IO::read_STL(filename_or_buf, mesh)) {
    return std::nullopt;
  }
  return mesh;
}

std::vector<Surface_mesh> SplitConnectedComponents(Surface_mesh& mesh) {
  using face_descriptor = boost::graph_traits<Surface_mesh>::face_descriptor;
  mesh.add_property_map<face_descriptor, Surface_mesh::faces_size_type>(
              "f:connected", Surface_mesh::faces_size_type(0));
  
  std::vector<Surface_mesh> components;
  PMP::split_connected_components(mesh, components);
  return components;
}

Plane_3 FindMiddleSlicePlane(const Surface_mesh& mesh) {
  const CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  const double mid_z = bbox.zmin() + bbox.z_span() / 2.;
  return Plane_3{Point_3{0, 0, mid_z}, Vector_3{0, 0, 1}};
}

Polylines SliceAtPlane(const Surface_mesh& mesh, Plane_3 plane) {
  using Polyline_3 = std::vector<Kernel::Point_3>;
  using Polylines_3 = std::vector<Polyline_3>;
  Polylines_3 polys_3;
  CGAL::Polygon_mesh_slicer<Surface_mesh, Kernel> slicer(mesh);
  slicer(plane, std::back_inserter(polys_3));
  Polylines polys_2;
  polys_2.reserve(polys_3.size());
  for (auto& poly_3 : polys_3) {
    Polyline poly_2;
    for (auto& p3 : poly_3) {
      Point_2 p2 = plane.to_2d(p3);
      poly_2.emplace_back(p2.x(), p2.y());
    }
    polys_2.push_back(poly_2);
  }
  return polys_2;
}

std::vector<std::tuple<double, double>> ComputeFacesOnZPlanes(
    const Surface_mesh& mesh, double tolerance = 1e-6) {
  std::vector<std::tuple<double, double>> result;
  std::unordered_map<double, double> z_plane_areas;

  for (Face_index f : mesh.faces()) {
    Vector_3 normal = PMP::compute_face_normal(f, mesh);
    if (std::fabs(normal.x()) <= tolerance &&
        std::fabs(normal.y()) <= tolerance) {
      auto vertices = CGAL::vertices_around_face(mesh.halfedge(f), mesh);
      double z_avg = 0.0;
      int count = 0;
      for (auto v : vertices) {
        z_avg += mesh.point(v).z();
        ++count;
      }
      z_avg /= count;
      double z_plane = std::round(z_avg / tolerance) * tolerance;
      double face_area = PMP::face_area(f, mesh);
      z_plane_areas[z_plane] += face_area;
    }
  }

  for (const auto& [z_plane, area] : z_plane_areas) {
    result.emplace_back(z_plane, area);
  }

  return result;
}

}  // namespace slice

namespace clean {

using namespace ClipperLib;
constexpr uint64_t kPrecisionMul = 100000l;
constexpr cInt toInt(double e) { return static_cast<int>(e * kPrecisionMul); }
constexpr double toDouble(cInt e) {
  return static_cast<double>(e) / kPrecisionMul;
}

ClipperLib::IntRect GetBounds(const Paths& paths) {
  ClipperLib::Clipper c_bounds;
  c_bounds.AddPaths(paths, ptSubject, true);
  c_bounds.StrictlySimple(true);
  Paths out;
  c_bounds.Execute(ctUnion, out, pftEvenOdd);
  return c_bounds.GetBounds();
}

Path PathFromDoublePoints(const std::vector<DoublePoint>& points) {
  Path path;
  path.reserve(points.size());
  for (const auto& item : points)
    path.emplace_back(toInt(item.X), toInt(item.Y));
  return path;
}

std::vector<DoublePoint> DoublePointsFromPath(const Path& path) {
  std::vector<DoublePoint> points;
  for (const auto& item : path)
    points.emplace_back(toDouble(item.X), toDouble(item.Y));
  return points;
}

Polylines PathsToPolylines(const Paths& paths) {
  Polylines poly;
  for (const auto& p : paths) {
    poly.push_back(DoublePointsFromPath(p));
  }
  return poly;
}

Paths PolylinesToPaths(const Polylines& polylines) {
  Paths paths;
  for (const auto& poly : polylines)
    paths.push_back(PathFromDoublePoints(poly));
  return paths;
}

ViewBox CleanPaths(Polylines& polylines) {
  Paths paths = PolylinesToPaths(polylines);
  Clipper c;
  c.AddPaths(paths, ptSubject, true);
  c.StrictlySimple(true);
  const IntRect bounds = c.GetBounds();
  c.Execute(ctUnion, paths, pftEvenOdd);
  polylines.clear();
  for (const auto& path : paths)
    polylines.push_back(DoublePointsFromPath(path));
  return {
      .minX = toDouble(bounds.left),
      .minY = toDouble(bounds.top),
      .width = toDouble(bounds.right - bounds.left),
      .height = toDouble(bounds.bottom - bounds.top),
  };
}

ViewBox ViewBoxOfPaths(const Paths& paths) {
  Clipper c;
  c.AddPaths(paths, ptSubject, true);
  const IntRect bounds = c.GetBounds();
  ViewBox vbox{
      .minX = clean::toDouble(bounds.left),
      .minY = clean::toDouble(bounds.top),
      .width = clean::toDouble(bounds.right - bounds.left),
      .height = clean::toDouble(bounds.bottom - bounds.top),
  };
  return vbox;
}

Paths Union(const Paths& paths) {
  Paths out;
  Clipper c;
  c.AddPaths(paths, ptSubject, true);
  c.StrictlySimple(true);
  c.Execute(ctUnion, out, pftNonZero);
  return out;
}

}  // namespace clean

namespace svg {

struct PolyDepth {
  Polylines line;
  double depth;
};
using PolyDepths = std::vector<PolyDepth>;

void WriteSvgD(const Polylines& polygon, std::stringstream& s) {
  auto Point = [&s](const ClipperLib::DoublePoint& pt) {
    s << pt.X << "," << pt.Y << "";
  };
  for (const auto& polyline : polygon) {
    auto it = polyline.cbegin();
    s << "M";
    Point(*it);
    ++it;
    for (; it != polyline.cend(); ++it) {
      s << " L";
      Point(*it);
    }
    s << " z ";
  }
}

void WritePath(const Polylines& polygon, const ViewBox& vbox,
               std::stringstream& s) {
  s << R"(<path d=")";
  WriteSvgD(polygon, s);
  s << R"(")"
    << R"( fill-rule="evenodd" stroke="black" stroke-width="1" fill="#333")"
    << R"( fill-opacity="0.25" vector-effect="non-scaling-stroke" shape-rendering="crispEdges"></path>)";
}

void WritePath(const PolyDepth& polygon, const ViewBox& vbox,
               std::stringstream& s) {
  s << R"(<path d=")";
  WriteSvgD(polygon.line, s);
  s << R"("></path>)";
}

std::string WriteSvgSplit(const PolyDepths& polygons, const ViewBox vbox) {
  std::stringstream s;
  s << R"(<svg xmlns="http://www.w3.org/2000/svg" viewBox=")"  //
    << vbox.minX << " " << vbox.minY << " " << vbox.width << " " << vbox.height
    << R"(" width=")" << vbox.width << R"(mm" height=")" << vbox.height
    << R"(mm")"
    << R"(><title property="dc:title">stltosvg v1.0 -- copyright (c) Alexandre Macabies</title><desc property="dc:creator">stltosvg v1.0 -- copyright (c) Alexandre Macabies</desc>)";
  for (const auto& poly : polygons) {
    WritePath(poly, vbox, s);
  }
  s << "</svg>";
  return s.str();
}

std::string WriteSvgEasel(const PolyDepths& polygons, const ViewBox vbox) {
  std::stringstream s;
  s << R"(<svg xmlns="http://www.w3.org/2000/svg" viewBox=")"  //
    << vbox.minX << " " << vbox.minY << " " << vbox.width << " " << vbox.height
    << R"(" width=")" << vbox.width << R"(mm" height=")" << vbox.height
    << R"(mm")"
    << R"(><title property="dc:title">stltosvg v1.0 -- copyright (c) Alexandre Macabies</title><desc property="dc:creator">stltosvg v1.0 -- copyright (c) Alexandre Macabies</desc>)"
    << R"(<g><g>)"
    // Bounding Box
    << R"s(<g stroke="none" fill="rgb(255,255,255)"><g><path d=")s"
    << "M" << vbox.minX << "," << vbox.minY                                  //
    << " L" << (vbox.minX + vbox.width) << "," << vbox.minY                  //
    << " L" << (vbox.minX + vbox.width) << "," << (vbox.minY + vbox.height)  //
    << " L" << (vbox.minX) << "," << (vbox.minY + vbox.height)               //
    << " z"
    << R"("></path></g></g>)";
  for (const auto& poly : polygons) {
    const uint8_t gray = static_cast<uint8_t>(255.0 * poly.depth);
    const int g = gray;
    s << R"s(<g stroke="none" fill="rgb()s" << g << "," << g << "," << g
      << R"s()"><g>)s";
    WritePath(poly, vbox, s);
    s << R"("></g></g>)";
  }
  s << R"(</g></g></svg>)";
  return s.str();
}

std::string WriteSvg(const Polylines& polygons, const ViewBox vbox) {
  std::stringstream s;
  s << R"(<svg xmlns="http://www.w3.org/2000/svg" viewBox=")"  //
    << vbox.minX << " " << vbox.minY << " " << vbox.width << " " << vbox.height
    << R"(" width=")" << vbox.width << R"(mm" height=")" << vbox.height
    << R"(mm")"
    << R"(><title property="dc:title">stltosvg v1.0 -- copyright (c) Alexandre Macabies</title><desc property="dc:creator">stltosvg v1.0 -- copyright (c) Alexandre Macabies</desc>)";
  WritePath(polygons, vbox, s);
  s << "</svg>";
  return s.str();
}

}  // namespace svg

// Logic moved outside #ifdef __EMSCRIPTEN__ for reusability
std::vector<std::string> StlToPaths(const std::string& stl, bool reorient) {
  std::istringstream is(stl);
  auto maybe_mesh = slice::ReadSTL(is);
  if (!maybe_mesh.has_value()) {
    return {};
  }

  auto mesh = reorient ? OrientModel(*maybe_mesh) : *maybe_mesh;

  std::vector<std::string> out;
  const auto comps = slice::SplitConnectedComponents(mesh);
  LOG("Found {} components", comps.size());
  for (const auto& cmesh : comps) {
    const auto plane = slice::FindMiddleSlicePlane(cmesh);
    Polylines polylines = slice::SliceAtPlane(cmesh, plane);
    ViewBox vbox = clean::CleanPaths(polylines);
    out.push_back(svg::WriteSvg(polylines, vbox));
  }
  return out;
}

std::string StlToEaselSvg(const std::string& stl, double area_tol, double nudge,
                          bool reorient, bool reverseOrder, bool reverseDepth) {
  std::istringstream is(stl);
  auto maybe_mesh = slice::ReadSTL(is);
  if (!maybe_mesh.has_value()) {
    return "";
  }

  auto mesh = reorient ? OrientModel(*maybe_mesh) : *maybe_mesh;

  const CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  LOG("Height: {} -> {} = {}mm", bbox.zmin(), bbox.zmax(), bbox.z_span());

  auto faces_on_z_planes = slice::ComputeFacesOnZPlanes(mesh);
  LOG("Found {} faces on z-planes", faces_on_z_planes.size());

  std::stable_sort(faces_on_z_planes.begin(), faces_on_z_planes.end(),
                   [reverseOrder](const std::tuple<double, double>& a,
                                  const std::tuple<double, double>& b) {
                     bool lower = std::get<0>(a) < std::get<0>(b);
                     return reverseOrder ? !lower : lower;
                   });

  ClipperLib::Paths bound_paths;
  std::vector<std::tuple<double, ClipperLib::Paths>> out_paths;
  // First pass: slice & collect bounds.
  for (auto [z, area] : faces_on_z_planes) {
    if (area <= 100 || IsEqual(z, 0.0, 1e-4)) continue;
    Polylines polylines = slice::SliceAtPlane(
        mesh, slice::Plane_3{slice::Point_3{0, 0, z + nudge}, slice::Vector_3{0, 0, 1}});
    LOG("Slicing at z={} (nudge={})", z, nudge);
    ClipperLib::Paths paths = clean::PolylinesToPaths(polylines);
    bound_paths.insert(bound_paths.end(), paths.cbegin(), paths.cend());
    out_paths.push_back({z, paths});
  }
  const auto viewbox = clean::ViewBoxOfPaths(bound_paths);
  const auto boundArea =
      clean::toInt(viewbox.height) * clean::toInt(viewbox.width);
  // Second pass: remove large areas meant for origin-alignment.
  for (auto&& [z, paths] : out_paths) {
    for (auto it = paths.begin(); it != paths.end();) {
      const auto area = std::abs(ClipperLib::Area(*it));
      if (area_tol != 1.0 && (area >= (double)boundArea * area_tol)) {
        it = paths.erase(it);
      } else if (area <= clean::toInt(1)) {
        it = paths.erase(it);
      } else {
        ++it;
      }
    }
  }
  // Third pass: nonzero (union) each layer.
  svg::PolyDepths depths;
  for (auto [z, paths] : out_paths) {
    if (paths.size() == 0) continue;
    double scaled_z = (z - bbox.zmin()) / bbox.z_span();
    double depth = reverseDepth ? scaled_z : 1.0 - scaled_z;
    depths.push_back({
        .line = clean::PathsToPolylines(clean::Union(paths)),
        .depth = depth,
    });
  }
  return svg::WriteSvgEasel(depths, viewbox);
}

#ifdef __EMSCRIPTEN__
using namespace emscripten;
EMSCRIPTEN_BINDINGS(Module) {
  register_vector<std::string>("Paths");
  function("StlToPaths", &StlToPaths);
  function("StlToEaselSvg", &StlToEaselSvg);
}
#endif