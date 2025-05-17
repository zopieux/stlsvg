#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>  // To ensure meshes are triangulated
#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/iterator.h>            // For face_iterator
#include <CGAL/linear_least_squares_fitting_2.h>  // For PCA-like behavior
#include <CGAL/number_utils.h>  // For CGAL::scalar_product, CGAL::cross_product, and potentially PI related constants

#include <algorithm>  // For std::max, std::min
#include <cmath>  // For std::sqrt, std::acos, std::abs, M_PI (if _USE_MATH_DEFINES)
#include <iostream>
#include <limits>  // For std::numeric_limits
#include <set>
#include <vector>

namespace {

// Define Kernel and Surface_mesh
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = Kernel::Vector_3;
using Vector_2 = Kernel::Vector_2;
using Point_2 = Kernel::Point_2;
using Line_2 = Kernel::Line_2;
using Direction_2 = Kernel::Direction_2;
using Plane_3 = Kernel::Plane_3;
using Bbox_3 = CGAL::Bbox_3;

using Polygon_2 = CGAL::Polygon_2<Kernel>;
using PolygonWithHoles_2 = CGAL::Polygon_with_holes_2<Kernel>;
using PolygonSet_2 = CGAL::Polygon_set_2<Kernel>;
using Polyline_2 = std::vector<Point_2>;
using Polyline_3 = std::vector<Point_3>;
using Polylines_3 = std::vector<Polyline_3>;

using Surface_mesh = CGAL::Surface_mesh<Point_3>;
using Aff_transformation_3 = CGAL::Aff_transformation_3<Kernel>;
using face_descriptor = Surface_mesh::Face_index;
using vertex_descriptor = Surface_mesh::Vertex_index;
namespace PMP = CGAL::Polygon_mesh_processing;

// Helper to normalize a vector, returns false if zero length
bool normalize_vector(Vector_3& v, double tolerance = 1e-9) {
  double sq_len = v.squared_length();
  if (sq_len < tolerance * tolerance) {  // Compare with squared tolerance
    v = Vector_3(0, 0, 0);               // Explicitly set to zero if too small
    return false;                        // Zero vector or too small
  }
  v = v / std::sqrt(sq_len);
  return true;
}

Aff_transformation_3 rotation_around_axis(const Vector_3& axis,
                                          double angle_rad) {
  typedef typename Kernel::RT RT;

  // Normalize the axis vector
  RT axis_length = CGAL::sqrt(axis.squared_length());
  if (axis_length == 0) {
    // Handle the case of a zero vector axis (no rotation)
    return CGAL::Aff_transformation_3<Kernel>(CGAL::IDENTITY);
  }
  CGAL::Vector_3<Kernel> normalized_axis = axis / axis_length;

  RT ux = normalized_axis.x();
  RT uy = normalized_axis.y();
  RT uz = normalized_axis.z();

  RT cos_theta = std::cos(angle_rad);
  RT sin_theta = std::sin(angle_rad);
  RT one_minus_cos_theta = 1.0 - cos_theta;

  // Construct the rotation matrix elements
  RT m00 = cos_theta + ux * ux * one_minus_cos_theta;
  RT m01 = ux * uy * one_minus_cos_theta - uz * sin_theta;
  RT m02 = ux * uz * one_minus_cos_theta + uy * sin_theta;
  RT m03 = 0.0;

  RT m10 = uy * ux * one_minus_cos_theta + uz * sin_theta;
  RT m11 = cos_theta + uy * uy * one_minus_cos_theta;
  RT m12 = uy * uz * one_minus_cos_theta - ux * sin_theta;
  RT m13 = 0.0;

  RT m20 = uz * ux * one_minus_cos_theta - uy * sin_theta;
  RT m21 = uz * uy * one_minus_cos_theta + ux * sin_theta;
  RT m22 = cos_theta + uz * uz * one_minus_cos_theta;
  RT m23 = 0.0;

  // Create the Aff_transformation_3 from the matrix
  return Aff_transformation_3(m00, m01, m02, m03, m10, m11, m12, m13, m20, m21,
                              m22, m23);
}

// Helper to get rotation to align 'from_vec' to 'to_vec'
Aff_transformation_3 get_rotation_align_vectors(Vector_3 from_vec,
                                                Vector_3 to_vec,
                                                double tolerance = 1e-9) {
  if (!normalize_vector(from_vec, tolerance) ||
      !normalize_vector(to_vec, tolerance)) {
    return Aff_transformation_3(CGAL::IDENTITY);  // Cannot align zero vector
  }

  double cos_angle = CGAL::scalar_product(from_vec, to_vec);
  // Clamp cos_angle to [-1, 1] to avoid domain errors with acos due to floating
  // point issues
  cos_angle = std::max(-1.0, std::min(1.0, cos_angle));

  // Check for collinear vectors
  if (std::abs(cos_angle) >=
      1.0 - tolerance * tolerance) {  // Using a tolerance related to cos_angle
                                      // for collinearity
    if (cos_angle > 0) {              // Already aligned (or very close)
      return Aff_transformation_3(CGAL::IDENTITY);
    } else {  // Opposite direction, 180-degree rotation
      // Find an arbitrary axis perpendicular to from_vec
      Vector_3 rot_axis = CGAL::cross_product(from_vec, Vector_3(1, 0, 0));
      if (rot_axis.squared_length() <
          tolerance * tolerance) {  // from_vec was parallel to X-axis
        rot_axis = CGAL::cross_product(from_vec, Vector_3(0, 1, 0));
      }
      // If from_vec is also parallel to Y (so it's Z-aligned), cross with
      // (0,1,0) is fine. Normalization for rot_axis (it must be unit for the
      // constructor)
      if (!normalize_vector(rot_axis, tolerance)) {
        // This should ideally not happen if from_vec is non-zero.
        // If from_vec is (0,0,1), cross with (1,0,0) -> (0,-1,0). Normalized.
        // OK. If from_vec is (1,0,0), cross with (1,0,0) -> (0,0,0). Then cross
        // with (0,1,0) -> (0,0,1). Normalized. OK. As a last resort if somehow
        // rot_axis is still zero (e.g. from_vec was (0,0,0) despite prior
        // normalization): Fallback: try a default axis like Z, or X if from_vec
        // is Z
        if (from_vec == Vector_3(0, 0, 1) || from_vec == Vector_3(0, 0, -1))
          rot_axis = Vector_3(1, 0, 0);
        else
          rot_axis = Vector_3(0, 0, 1);
        // No need to normalize these canonical axes again
      }
      return rotation_around_axis(rot_axis, CGAL_PI);
      // return Aff_transformation_3(CGAL::ROTATION, rot_axis, 180.0); // Angle
      // in degrees
    }
  }

  // General case: vectors are not collinear
  Vector_3 rot_axis = CGAL::cross_product(from_vec, to_vec);
  // The cross product vector's length is |from_vec||to_vec|sin(angle).
  // Since from_vec and to_vec are unit, length is sin(angle).
  // rot_axis itself must be normalized for the constructor.
  if (!normalize_vector(rot_axis, tolerance)) {
    // This case should not be reached if vectors are not collinear and
    // non-zero. If it is, means something is wrong, possibly numerical
    // instability.
    return Aff_transformation_3(
        CGAL::IDENTITY);  // Fallback, should be investigated if hit
  }

  double angle_rad = std::acos(cos_angle);  // Angle is in [0, PI]
  return rotation_around_axis(rot_axis, angle_rad);
}

// Custom comparator for Vector_3 for std::set (based on rounded values for
// uniqueness)
struct Vector3ApproxComparator {
  double tol_sq;  // Use squared tolerance for squared length comparisons if
                  // needed, direct for components
  double tol;
  Vector3ApproxComparator(double t = 1e-5) : tol(t), tol_sq(t * t) {}
  bool operator()(const Vector_3& a, const Vector_3& b) const {
    // A more robust way for "equality" check before comparison for ordering:
    if (CGAL::squared_distance(Point_3(a.x(), a.y(), a.z()),
                               Point_3(b.x(), b.y(), b.z())) < tol_sq) {
      return false;  // Treat as equal for set uniqueness if very close
    }
    // Lexicographical comparison
    if (std::abs(a.x() - b.x()) > tol) return a.x() < b.x();
    if (std::abs(a.y() - b.y()) > tol) return a.y() < b.y();
    if (std::abs(a.z() - b.z()) > tol) return a.z() < b.z();
    return false;  // Fallback: consider equal if all components are within tol
  }
};

}

Surface_mesh OrientModel(Surface_mesh mesh, double tolerance = 1e-4) {
  if (CGAL::is_empty(mesh) || num_faces(mesh) == 0) {
    return mesh;  // Nothing to do
  }

  std::vector<Vector_3> all_face_normals;
  all_face_normals.reserve(num_faces(mesh));

  Surface_mesh::Property_map<face_descriptor, Vector_3> fnormals_map;
  bool created;
  // Check if property map already exists, otherwise add it.
  auto opt_fnormals_map =
      mesh.property_map<face_descriptor, Vector_3>("f:normal");
  if (opt_fnormals_map.second) {  // .second is true if map exists
    fnormals_map = opt_fnormals_map.first;
  } else {
    boost::tie(fnormals_map, created) =
        mesh.add_property_map<face_descriptor, Vector_3>("f:normal",
                                                         Vector_3(0, 0, 0));
  }
  PMP::compute_face_normals(mesh, fnormals_map,
                            PMP::parameters::geom_traits(Kernel()));

  for (face_descriptor fd : faces(mesh)) {
    Vector_3 n = fnormals_map[fd];  // Already computed normal
    if (normalize_vector(
            n, tolerance)) {  // Normalize and check if it's not a zero vector
      all_face_normals.push_back(n);
    }
  }

  if (all_face_normals.empty()) {
    std::cerr << "Warning: No valid face normals found." << std::endl;
    return mesh;  // Or apply identity if no normals to guide
  }

  std::set<Vector_3, Vector3ApproxComparator> unique_model_normals_set{
      Vector3ApproxComparator(tolerance)};
  for (const auto& n : all_face_normals) {
    unique_model_normals_set.insert(n);
    unique_model_normals_set.insert(
        -n);  // Also consider the opposite direction
  }

  std::vector<Vector_3> candidate_model_up_vectors(
      unique_model_normals_set.begin(), unique_model_normals_set.end());
  if (candidate_model_up_vectors.empty()) {
    candidate_model_up_vectors.push_back(
        Vector_3(0, 0, 1));  // Default fallback
  }

  Aff_transformation_3 best_transform(CGAL::IDENTITY);
  long best_score_undercuts =
      std::numeric_limits<long>::max();  // Want to minimize
  long best_score_bottoms = -1;          // Want to maximize
  long best_score_walls = -1;            // Want to maximize

  Vector_3 cnc_z_axis(0, 0, 1);
  Vector_3 cnc_x_axis(1, 0, 0);

  for (const Vector_3& model_up_candidate : candidate_model_up_vectors) {
    Aff_transformation_3 primary_rotation =
        get_rotation_align_vectors(model_up_candidate, cnc_z_axis, tolerance);

    std::vector<Vector_3> model_secondary_candidates_bases;
    Vector_3 model_x_trial_world_frame = CGAL::cross_product(
        cnc_z_axis,
        Vector_3(0, 1, 0));  // Try to align model's "something" to world X
                             // This is world Y if cnc_z is Z.
                             // This part needs to be model relative.

    // Candidates for model's local X axis (must be perpendicular to
    // model_up_candidate)
    Vector_3 model_local_x_candidate =
        CGAL::cross_product(model_up_candidate, Vector_3(0, 1, 0));
    if (!normalize_vector(model_local_x_candidate, tolerance)) {
      model_local_x_candidate =
          CGAL::cross_product(model_up_candidate, Vector_3(1, 0, 0));
      normalize_vector(model_local_x_candidate,
                       tolerance);  // Ensure it's valid
    }

    if (model_local_x_candidate.squared_length() > tolerance * tolerance) {
      model_secondary_candidates_bases.push_back(model_local_x_candidate);
    }
    // Also try canonical model axes as candidates for what becomes CNC X after
    // primary rotation
    model_secondary_candidates_bases.push_back(
        primary_rotation(Vector_3(1, 0, 0)));
    model_secondary_candidates_bases.push_back(
        primary_rotation(Vector_3(0, 1, 0)));
    // Filter out any candidates that are too close to cnc_z_axis after primary
    // rotation
    model_secondary_candidates_bases.erase(
        std::remove_if(
            model_secondary_candidates_bases.begin(),
            model_secondary_candidates_bases.end(),
            [&](const Vector_3& v_transformed) {
              Vector_3 v_proj_xy(v_transformed.x(), v_transformed.y(), 0);
              return v_proj_xy.squared_length() <
                     tolerance *
                         tolerance;  // Remove if projection is near zero
            }),
        model_secondary_candidates_bases.end());

    if (model_secondary_candidates_bases
            .empty()) {  // If all candidates were bad or Z-aligned
      model_secondary_candidates_bases.push_back(
          Vector_3(1, 0, 0));  // Default to world X
    }

    for (const Vector_3& model_axis_to_align_with_cnc_x_after_primary_rot :
         model_secondary_candidates_bases) {
      // This vector is already in the frame *after* primary_rotation
      // We want to align its XY projection with cnc_x_axis
      Vector_3 vec_proj_cnc_xy(
          model_axis_to_align_with_cnc_x_after_primary_rot.x(),
          model_axis_to_align_with_cnc_x_after_primary_rot.y(), 0.0);

      Aff_transformation_3 secondary_rotation(CGAL::IDENTITY);
      if (normalize_vector(vec_proj_cnc_xy,
                           tolerance)) {  // only if projection is non-zero
        secondary_rotation =
            get_rotation_align_vectors(vec_proj_cnc_xy, cnc_x_axis, tolerance);
      }

      Aff_transformation_3 current_transform =
          secondary_rotation * primary_rotation;

      long current_undercuts = 0;
      long current_bottoms = 0;
      long current_walls = 0;

      for (const Vector_3& original_normal : all_face_normals) {
        Vector_3 n_world = current_transform(original_normal);
        if (!normalize_vector(n_world, tolerance)) continue;

        double z_comp = n_world.z();
        // Use a slightly larger tolerance for "is vertical/horizontal/up/down"
        // than for vector normalization or geometric calculations.
        double align_tol = 1e-3;  // Tolerance for alignment checks

        if (z_comp > 1.0 - align_tol) {  // Normal is [0,0,1] (pointing up)
          current_bottoms++;
        } else if (z_comp <
                   -1.0 + align_tol) {  // Normal is [0,0,-1] (pointing down)
          current_undercuts++;
        } else if (std::abs(z_comp) <
                   align_tol) {  // Normal is in XY plane (vertical wall)
          current_walls++;
        } else if (z_comp < -align_tol) {  // General undercut (pointing
                                           // somewhat downwards)
          current_undercuts++;
        }
        // else: normal is pointing somewhat upwards, not a bottom, not a wall,
        // not an undercut.
      }

      bool update_best = false;
      if (current_undercuts < best_score_undercuts) {
        update_best = true;
      } else if (current_undercuts == best_score_undercuts) {
        if (current_bottoms > best_score_bottoms) {
          update_best = true;
        } else if (current_bottoms == best_score_bottoms) {
          if (current_walls > best_score_walls) {
            update_best = true;
          }
        }
      }

      if (update_best) {
        best_transform = current_transform;
        best_score_undercuts = current_undercuts;
        best_score_bottoms = current_bottoms;
        best_score_walls = current_walls;
      }
    }
  }

  PMP::transform(best_transform, mesh);

  return mesh;
}

double get_xy_projected_aabb_area(const Surface_mesh& mesh) {
  if (CGAL::is_empty(mesh) || num_vertices(mesh) == 0) {
    return std::numeric_limits<double>::max();  // Return large value for
                                                // empty/invalid meshes
  }

  auto vpmap = get(CGAL::vertex_point, mesh);  // Property map for vertex points
  double min_x = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double min_y = std::numeric_limits<double>::max();
  double max_y = std::numeric_limits<double>::lowest();

  bool has_points = false;
  for (Surface_mesh::Vertex_index vi : vertices(mesh)) {
    Point_3 p = vpmap[vi];
    // Use CGAL::to_double for robust conversion from Kernel::FT to double
    min_x = std::min(min_x, CGAL::to_double(p.x()));
    max_x = std::max(max_x, CGAL::to_double(p.x()));
    min_y = std::min(min_y, CGAL::to_double(p.y()));
    max_y = std::max(max_y, CGAL::to_double(p.y()));
    has_points = true;
  }

  if (!has_points) return std::numeric_limits<double>::max();
  // Handle cases where all points are collinear (resulting in zero width or
  // height)
  if (max_x < min_x || max_y < min_y) return 0.0;

  return (max_x - min_x) * (max_y - min_y);
}

Surface_mesh AlignMeshToXYAxes(
    Surface_mesh mesh) {  // Input mesh is by value (a copy)
  if (CGAL::is_empty(mesh) || num_vertices(mesh) < 2) {
    // Need at least 2 distinct points for PCA to define a line.
    // For <3 non-collinear points, 2D AABB might not be unique.
    return mesh;
  }

  // 1. Project Vertices to XY Plane
  std::vector<Point_2> projected_points;
  projected_points.reserve(num_vertices(mesh));
  auto vpmap = get(CGAL::vertex_point, mesh);
  for (Surface_mesh::Vertex_index vi : vertices(mesh)) {
    Point_3 p3 = vpmap[vi];
    projected_points.emplace_back(p3.x(), p3.y());
  }

  // Ensure we have enough distinct points for fitting
  // (CGAL::linear_least_squares_fitting_2 might handle this, but good to check)
  if (projected_points.size() < 2) {
    return mesh;
  }
  // Check for distinct points (optional, fitting might handle degenerate cases)
  std::sort(projected_points.begin(), projected_points.end());
  projected_points.erase(
      std::unique(projected_points.begin(), projected_points.end()),
      projected_points.end());
  if (projected_points.size() < 2) {
    return mesh;  // Not enough distinct points
  }

  // 2. Compute Principal Component Axis
  Line_2 principal_axis_line;
  // The last argument CGAL::Dimension_tag<0>() means we are fitting points.
  CGAL::linear_least_squares_fitting_2(
      projected_points.begin(), projected_points.end(), principal_axis_line,
      CGAL::Dimension_tag<0>());

  Direction_2 pca_direction_cgal = principal_axis_line.direction();
  Vector_2 pca_dir_2d(pca_direction_cgal.dx(), pca_direction_cgal.dy());

  // If PCA direction is near zero (e.g., all points very close or collinear in
  // a way fitting fails)
  if (pca_dir_2d.squared_length() < 1e-12) {
    return mesh;  // Cannot determine a principal direction robustly
  }

  // 3. Candidate Rotations
  // Angle of pca_dir_2d relative to the positive X-axis
  double pca_angle_rad = std::atan2(CGAL::to_double(pca_dir_2d.y()),
                                    CGAL::to_double(pca_dir_2d.x()));

  std::vector<double> candidate_rotation_angles_rad;
  // Angle to rotate the *mesh* so that pca_dir aligns with X
  candidate_rotation_angles_rad.push_back(-pca_angle_rad);
  // Angle to rotate the *mesh* so that pca_dir aligns with Y (or its
  // perpendicular aligns with X)
  candidate_rotation_angles_rad.push_back(-pca_angle_rad + M_PI / 2.0);
  // No need to test +PI or -PI/2 versions due to AABB area symmetry for 180-deg
  // flips.

  double best_final_rotation_angle_rad =
      0.0;  // Assume initially no rotation is best
  double min_aabb_area =
      get_xy_projected_aabb_area(mesh);  // Area of the original orientation

  // 4. Evaluate and Select
  for (double rot_angle_rad : candidate_rotation_angles_rad) {
    // Normalize the test angle to avoid very large angle values if accumulated
    double current_test_angle_rad = std::fmod(rot_angle_rad, 2.0 * M_PI);
    if (current_test_angle_rad > M_PI) current_test_angle_rad -= 2.0 * M_PI;
    if (current_test_angle_rad <= -M_PI) current_test_angle_rad += 2.0 * M_PI;

    // If this normalized angle is very close to zero, skip (already covered by
    // initial state)
    if (std::abs(current_test_angle_rad) < 1e-7 &&
        best_final_rotation_angle_rad == 0.0) {
      continue;
    }

    Surface_mesh temp_mesh =
        mesh;  // Start with a fresh copy of the *original input* Z-up mesh

    Aff_transformation_3 z_rotation =
        rotation_around_axis(Vector_3(0, 0, 1), current_test_angle_rad);

    PMP::transform(z_rotation, temp_mesh);
    double current_aabb_area = get_xy_projected_aabb_area(temp_mesh);

    // Use a small tolerance (epsilon) for comparing areas to avoid floating
    // point noise
    if (current_aabb_area < min_aabb_area - 1e-9) {
      min_aabb_area = current_aabb_area;
      best_final_rotation_angle_rad =
          current_test_angle_rad;  // Store the effective angle
    }
  }

  // 5. Apply the best rotation to the input mesh (which is 'mesh' here, being a
  // copy)
  if (std::abs(best_final_rotation_angle_rad) >
      1e-7) {  // Apply if a meaningful rotation was found
    Aff_transformation_3 final_z_rotation =
        rotation_around_axis(Vector_3(0, 0, 1), best_final_rotation_angle_rad);
    PMP::transform(final_z_rotation, mesh);
  }

  return mesh;
}

#if 0
int main(int argc, char* argv[]) {
  const std::string filename =
      (argc > 1) ? argv[1] : "input.stl";  // Default to input.stl
  Surface_mesh mesh;
  std::ifstream stl_stream(filename, std::ios::in | std::ios::binary);
  if (!stl_stream) {
    std::cerr << "Error: Cannot open file: " << filename << std::endl;
    return 1;
  }
  if (!CGAL::IO::read_STL(stl_stream, mesh) || CGAL::is_empty(mesh) ||
      !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Error: Cannot read STL file or mesh is not valid/empty: "
              << filename << std::endl;
    return 1;
  }
  std::cout << "Input mesh: " << num_vertices(mesh) << " vertices, "
            << num_faces(mesh) << " faces." << std::endl;

  Surface_mesh oriented_mesh = OrientModel(mesh, 1e-4);  // Pass tolerance

  std::cout << "Oriented mesh: " << num_vertices(oriented_mesh) << " vertices, "
            << num_faces(oriented_mesh) << " faces." << std::endl;

  Surface_mesh rotated_mesh = AlignMeshToXYAxes(oriented_mesh);

  std::string output_filename = "/tmp/oriented_output.stl";
  if (!CGAL::IO::write_STL(output_filename, oriented_mesh,
                           CGAL::parameters::stream_precision(17))) {
    std::cerr << "Error writing STL file: " << output_filename << std::endl;
    return 1;
  }
  std::cout << "Oriented mesh written to " << output_filename << std::endl;

  output_filename = "/tmp/oriented_output_rotated.stl";
  if (!CGAL::IO::write_STL(output_filename, rotated_mesh,
                           CGAL::parameters::stream_precision(17))) {
    std::cerr << "Error writing STL file: " << output_filename << std::endl;
    return 1;
  }
  std::cout << "Rotated mesh written to " << output_filename << std::endl;

  return 0;
}
#endif
