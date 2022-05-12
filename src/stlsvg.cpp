#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Surface_mesh.h>

#ifdef __EMSCRIPTEN__
#include <sstream>
#endif

#include "clipper/clipper.hpp"

namespace {

using Polyline = std::vector<ClipperLib::DoublePoint>;
using Polylines = std::vector<Polyline>;

struct ViewBox {
  double minX, minY, width, height;
};

}  // namespace

namespace slice {

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Plane_3 = Kernel::Plane_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
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

}  // namespace slice

namespace clean {

using namespace ClipperLib;
constexpr uint64_t kPrecisionMul = 100000l;
constexpr cInt toInt(double e) { return static_cast<int>(e * kPrecisionMul); }
constexpr double toDouble(cInt e) {
  return static_cast<double>(e) / kPrecisionMul;
}

Path PathFromDoublePoints(const std::vector<DoublePoint>& points) {
  Path path;
  path.reserve(points.size());
  for (const auto& item : points)
    path.emplace_back(toInt(item.X), toInt(item.Y));
  return path;
}

std::vector<DoublePoint> DoublePointsFromPath(const ClipperLib::Path& path) {
  std::vector<DoublePoint> points;
  for (const auto& item : path)
    points.emplace_back(toDouble(item.X), toDouble(item.Y));
  return points;
}

ViewBox CleanPaths(Polylines& polylines) {
  Paths paths;
  for (const auto& poly : polylines)
    paths.push_back(PathFromDoublePoints(poly));
  Clipper c;
  c.AddPaths(paths, ClipperLib::ptSubject, true);
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

}  // namespace clean

namespace svg {

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
    // << R"( data-viewbox=")" << vbox.minX << " " << vbox.minY << " " <<
    // vbox.width << " " << vbox.height
    // << R"(" data-width=")" << vbox.width << R"(mm" data-height=")" <<
    // vbox.height << R"(mm")"
    << R"( fill-rule="evenodd" stroke="black" stroke-width="1" fill="#333")"
    << R"( fill-opacity="0.25" vector-effect="non-scaling-stroke" shape-rendering="crispEdges"></path>)";
}

std::string WriteSvg(const Polylines& polygon, const ViewBox vbox) {
  std::stringstream s;
  s << R"(<svg xmlns="http://www.w3.org/2000/svg" viewBox=")"  //
    << vbox.minX << " " << vbox.minY << " " << vbox.width << " " << vbox.height
    << R"(" width=")" << vbox.width << R"(mm" height=")" << vbox.height
    << R"(mm")"
    << R"(><title property="dc:title">stltosvg v1.0 -- copyright (c) Alexandre Macabies</title><desc property="dc:creator">stltosvg v1.0 -- copyright (c) Alexandre Macabies</desc>)";
  WritePath(polygon, vbox, s);
  s << R"#("</svg>)#";
  return s.str();
}

}  // namespace svg

#ifdef __EMSCRIPTEN__
std::vector<std::string> StlToPaths(const std::string& stl) {
  std::vector<std::string> out;
  std::istringstream is(stl);
  auto maybe_mesh = slice::ReadSTL(is);
  if (!maybe_mesh.has_value()) {
    return out;
  }
  const auto comps = slice::SplitConnectedComponents(*maybe_mesh);
  int i = 0;
  // std::stringstream s;
  for (const auto& cmesh : comps) {
    const auto plane = slice::FindMiddleSlicePlane(cmesh);
    Polylines polylines = slice::SliceAtPlane(cmesh, plane);
    std::cerr << "component " << i++ << " has " << polylines.size()
              << " polylines" << std::endl;
    ViewBox vbox = clean::CleanPaths(polylines);
    // s.clear();
    out.push_back(svg::WriteSvg(polylines, vbox));
    // out.push_back(s.str());
  }
  return out;
}

#include <emscripten/bind.h>
using namespace emscripten;
EMSCRIPTEN_BINDINGS(Module) {
  register_vector<std::string>("Paths");
  // register_vector<std::vector<char>>("Buffers");
  function("StlToPaths", &StlToPaths);
}
#else
int main(int argc, char* argv[]) {
  const std::string filename = argv[1];
  auto opt_mesh = slice::ReadSTL(filename);
  if (!opt_mesh.has_value()) return 1;
  auto mesh = *opt_mesh;
  const auto comps = slice::SplitConnectedComponents(mesh);
  int i = 0;
  std::stringstream s;
  for (const auto& cmesh : comps) {
    const auto plane = slice::FindMiddleSlicePlane(cmesh);
    Polylines polylines = slice::SliceAtPlane(cmesh, plane);
    std::cerr << "component " << i++ << " has " << polylines.size()
              << " polylines" << std::endl;
    ViewBox viewBox = clean::CleanPaths(polylines);
    std::cout << svg::WriteSvg(polylines, viewBox) << "\n";
  }
  std::cout << s.str();
}
#endif
