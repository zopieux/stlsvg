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

#include "orient.cpp"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/bind.h>

#include <sstream>
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
  auto connected_map =
      mesh.add_property_map<face_descriptor, Surface_mesh::faces_size_type>(
              "f:connected", Surface_mesh::faces_size_type(0))
          .first;
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

#if 0
Point_3 compute_centroid(const Surface_mesh& mesh) {
    Vector_3 centroid(0, 0, 0);
    std::size_t vertex_count = 0;

    for (auto v : mesh.vertices()) {
        centroid = centroid + (mesh.point(v) - CGAL::ORIGIN);
        vertex_count++;
    }

    if (vertex_count > 0) {
        centroid = centroid / vertex_count;
    }

    return CGAL::ORIGIN + centroid;
}
// Function to scale a mesh by a given factor from its center
void scale_mesh(Surface_mesh& mesh, double scale_factor) {
    // Step 1: Compute the centroid of the mesh
    Point_3 centroid = compute_centroid(mesh);

    // Step 2: Construct the scaling transformation
    Aff_3 translate_to_origin(CGAL::TRANSLATION, -Vector_3(centroid - CGAL::ORIGIN));
    Aff_3 scaling(CGAL::SCALING, scale_factor);
    Aff_3 translate_back(CGAL::TRANSLATION, Vector_3(centroid - CGAL::ORIGIN));

    // Combined transformation: translate to origin -> scale -> translate back
    Aff_3 combined = translate_back * scaling * translate_to_origin;

    // Step 3: Apply the transformation to all points in the mesh
    PMP::transform(combined, mesh);
}

// Function to compute the bounding box of a mesh
Surface_mesh create_bounding_box(const Surface_mesh& mesh) {
  // Compute the bounding box
  CGAL::Bbox_3 bbox = PMP::bbox(mesh);

  // Create the bounding box corners
  Point_3 p_min(bbox.xmin(), bbox.ymin(), bbox.zmin());
  Point_3 p_max(bbox.xmax(), bbox.ymax(), bbox.zmax());

  // Generate a cuboid (bounding box) mesh
  Surface_mesh bbox_mesh;
  CGAL::make_hexahedron(p_min, Point_3(bbox.xmax(), bbox.ymin(), bbox.zmin()),
                        Point_3(bbox.xmax(), bbox.ymax(), bbox.zmin()),
                        Point_3(bbox.xmin(), bbox.ymax(), bbox.zmin()),
                        Point_3(bbox.xmin(), bbox.ymin(), bbox.zmax()),
                        Point_3(bbox.xmax(), bbox.ymin(), bbox.zmax()),
                        Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax()),
                        Point_3(bbox.xmin(), bbox.ymax(), bbox.zmax()),
                        bbox_mesh);

  return bbox_mesh;
}
#endif

std::vector<std::tuple<double, double>> ComputeFacesOnZPlanes(
    const Surface_mesh& mesh, double tolerance = 1e-6) {
  // Store the result as a vector of tuples
  std::vector<std::tuple<double, double>> result;

  // Map for grouping faces by z-coordinate planes
  // std::unordered_map<double, std::set<Face_index>> z_plane_faces;
  std::unordered_map<double, double> z_plane_areas;

  // Iterate over all faces in the mesh
  for (Face_index f : mesh.faces()) {
    // Compute the normal of the face
    Vector_3 normal = PMP::compute_face_normal(f, mesh);

    // Check if the normal is aligned with the z-axis
    if (std::fabs(normal.x()) <= tolerance &&
        std::fabs(normal.y()) <= tolerance) {
      // Get the face vertices
      auto vertices = CGAL::vertices_around_face(mesh.halfedge(f), mesh);

      // Compute the average z-coordinate of the face vertices
      double z_avg = 0.0;
      int count = 0;
      for (auto v : vertices) {
        z_avg += mesh.point(v).z();
        ++count;
      }
      z_avg /= count;

      // Round z_avg to the nearest plane within tolerance
      double z_plane = std::round(z_avg / tolerance) * tolerance;

      // Add the face to the appropriate z-plane group
      // z_plane_faces[z_plane].insert(f);

      // Compute the area of the face and add it to the total area for the plane
      double face_area = PMP::face_area(f, mesh);
      z_plane_areas[z_plane] += face_area;
    }
  }

  // Build the result from the maps
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
  s << R"#("</svg>)#";
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
    // Fucking BOUNDING BOX (depth=0 hack).
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
    s << R"("</g></g>)";
  }
  s << R"(</g></g>)" << R"#("</svg>)#";
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
  s << R"#("</svg>)#";
  return s.str();
}

}  // namespace svg

#ifdef __EMSCRIPTEN__
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
  int i = 0;
  for (const auto& cmesh : comps) {
    const auto plane = slice::FindMiddleSlicePlane(cmesh);
    Polylines polylines = slice::SliceAtPlane(cmesh, plane);
    ViewBox vbox = clean::CleanPaths(polylines);
    out.push_back(svg::WriteSvg(polylines, vbox));
  }
  return out;
}

#include <format>

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

  using clean::Path;
  using clean::Paths;

  Paths bound_paths;
  std::vector<std::tuple<double, Paths>> out_paths;
  // First pass: slice & collect bounds.
  for (auto [z, area] : faces_on_z_planes) {
    if (area <= 100 || IsEqual(z, 0.0, 1e-4)) continue;
    Polylines polylines = slice::SliceAtPlane(
        mesh, Plane_3{Point_3{0, 0, z + nudge}, Vector_3{0, 0, 1}});
    LOG("Slicing at z={} (nudge={})", z, nudge);
    Paths paths = clean::PolylinesToPaths(polylines);
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
      if (area_tol != 1.0 && (area >= boundArea * area_tol)) {
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

using namespace emscripten;
EMSCRIPTEN_BINDINGS(Module) {
  register_vector<std::string>("Paths");
  function("StlToPaths", &StlToPaths);
  function("StlToEaselSvg", &StlToEaselSvg);
}
#else

int main(int argc, char* argv[]) {
  const std::string filename = argv[1];
  auto opt_mesh = slice::ReadSTL(filename);
  if (!opt_mesh.has_value()) return 1;
  auto mesh = *opt_mesh;

  if (false) {
    // Define the transformation to swap Y and Z axes
    slice::Aff_3 swap_yz(  //
        1, 0, 0,           //
        0, 0, 1,           //
        0, 1, 0);
    for (auto pt = mesh.points().begin(); pt != mesh.points().end(); ++pt) {
      *pt = swap_yz(*pt);
    }
  }
  using namespace slice;
  using namespace clean;

  // ...
}
#endif
