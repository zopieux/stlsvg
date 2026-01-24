#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <getopt.h>
#include <filesystem>

// Forward declarations of functions from stlsvg.cpp
std::vector<std::string> StlToPaths(const std::string& stl, bool reorient);
std::string StlToEaselSvg(const std::string& stl, double area_tol, double nudge,
                          bool reorient, bool reverseOrder, bool reverseDepth);

void print_usage(const char* prog_name) {
    std::cerr << "Usage: " << prog_name << " <input.stl>\n"
              << "Options:\n"
              << "  --output <dir>       Output directory for SVG files\n"
              << "  --reorient           Reorient the model before slicing\n"
              << "  --reverse-order      Reverse visit order\n"
              << "  --reverse-depth      Reverse depth mapping\n"
              << "  --nudge <val>        Nudge value (default: -0.001)\n"
              << "  --help               Show this help message\n";
}

int main(int argc, char* argv[]) {
    std::string output_dir = "./test_output";
    bool reorient = false;
    bool reverse_order = false;
    bool reverse_depth = false;
    double nudge = -0.001;
    double area_tol = 1.0; // Default from stlsvg.js seems to be 1.0 passed to StlToEaselSvg

    struct option long_options[] = {
        {"output", required_argument, 0, 'o'},
        {"reorient", no_argument, 0, 'r'},
        {"reverse-order", no_argument, 0, 'O'},
        {"reverse-depth", no_argument, 0, 'D'},
        {"nudge", required_argument, 0, 'n'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "o:rODn:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'o':
                output_dir = optarg;
                break;
            case 'r':
                reorient = true;
                break;
            case 'O':
                reverse_order = true;
                break;
            case 'D':
                reverse_depth = true;
                break;
            case 'n':
                try {
                    nudge = std::stod(optarg);
                } catch (...) {
                    std::cerr << "Invalid nudge value: " << optarg << "\n";
                    return 1;
                }
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (optind >= argc) {
        std::cerr << "Error: Input file required.\n";
        print_usage(argv[0]);
        return 1;
    }

    std::string input_path = argv[optind];
    std::filesystem::path input_p(input_path);
    std::string stem = input_p.stem().string();

    // Read input file
    std::ifstream file(input_path, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open input file: " << input_path << "\n";
        return 1;
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string stl_content = buffer.str();

    // Ensure output directory exists
    try {
        std::filesystem::create_directories(output_dir);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating output directory: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Processing " << input_path << "...\n";

    // Generate Easel SVG
    std::string easel_svg = StlToEaselSvg(stl_content, area_tol, nudge, reorient, reverse_order, reverse_depth);
    if (!easel_svg.empty()) {
        std::filesystem::path easel_path = std::filesystem::path(output_dir) / (stem + "-easel.svg");
        std::ofstream out(easel_path);
        if (out) {
            out << easel_svg;
            std::cout << "Written " << easel_path << "\n";
        } else {
            std::cerr << "Error writing " << easel_path << "\n";
        }
    } else {
        std::cerr << "Failed to generate Easel SVG (possibly invalid STL)\n";
    }

    // Generate individual path SVGs
    std::vector<std::string> path_svgs = StlToPaths(stl_content, reorient);
    for (size_t i = 0; i < path_svgs.size(); ++i) {
        std::filesystem::path path_p = std::filesystem::path(output_dir) / (stem + "-subpart-" + std::to_string(i) + ".svg");
        std::ofstream out(path_p);
        if (out) {
            out << path_svgs[i];
            std::cout << "Written " << path_p << "\n";
        } else {
            std::cerr << "Error writing " << path_p << "\n";
        }
    }

    return 0;
}
