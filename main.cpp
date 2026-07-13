#include "db/PlacementDB.h"
#include "evaluator/HPWLEvaluator.h"
#include "parser/LimboBookshelfAdapter.h"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>



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
        HPWLEvaluator hpwl_evaluator;
        double total_hpwl = 0;
        for(std::size_t i = 0; i < db.nets().size(); i++){
            const int net_id = static_cast<int>(i);
            const double net_hpwl = hpwl_evaluator.netHPWL(db, net_id);
            db.addNetHPWL(net_id, net_hpwl);
            total_hpwl += net_hpwl;
        }

        std::ofstream fout("../result/placementdb_summary.txt");
        db.printSummary(fout);
        fout << "\n[Initial Placement Evaluation]\n";
        fout << std::fixed << std::setprecision(3);
        fout << "Total HPWL: " << total_hpwl << "\n";

        db.printSummary(std::cout);
        std::cout << "========== Initial Placement Evaluation ==========\n";
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "Total HPWL: " << total_hpwl << "\n";
        std::cout << "==================================================\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n'; return 2;
    }
    return 0;
}
