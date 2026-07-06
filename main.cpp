#include "db/PlacementDB.h"
#include "parser/LimboBookshelfAdapter.h"

#include <iostream>
#include <string>

namespace {
void printUsage(const char* argv0) {
    std::cerr << "Usage:\n  " << argv0 << " --aux <path_to_bookshelf_aux>\n";
}
}

int main(int argc, char** argv) {
    std::string aux_path;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--aux" && i + 1 < argc) aux_path = argv[++i];
        else if (arg == "--help" || arg == "-h") { printUsage(argv[0]); return 0; }
        else { std::cerr << "Unknown or incomplete argument: " << arg << '\n'; printUsage(argv[0]); return 1; }
    }
    if (aux_path.empty()) { printUsage(argv[0]); return 1; }
    try {
        PlacementDB db;
        LimboBookshelfAdapter adapter(db);
        if (!adapter.read(aux_path)) { std::cerr << "Failed to read Bookshelf design: " << aux_path << '\n'; return 2; }
        db.printSummary(std::cout);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n'; return 2;
    }
    return 0;
}
