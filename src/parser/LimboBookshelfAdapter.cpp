#include "parser/LimboBookshelfAdapter.h"
#include "db/PlacementDB.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(HAVE_LIMBO_BOOKSHELF)
#include <limbo/parsers/bookshelf/bison/BookshelfDriver.h>
namespace {
class CallbackDB : public BookshelfParser::BookshelfDataBase {
public:
    explicit CallbackDB(PlacementDB& db) : db_(db) {}
    void resize_bookshelf_node_terminals(int, int) override {}
    void resize_bookshelf_net(int) override {}
    void resize_bookshelf_pin(int) override {}
    void resize_bookshelf_row(int) override {}
    void add_bookshelf_terminal(std::string& name, int w, int h) override { db_.addCell(name, w, h, true); }
    void add_bookshelf_node(std::string& name, int w, int h) override { db_.addCell(name, w, h, false); }
    void add_bookshelf_net(BookshelfParser::Net const& net) override {
        const int net_id = db_.addNet(net.net_name);
        for (const auto& p : net.vNetPin) {
            if (!db_.hasCell(p.node_name)) throw std::runtime_error(".nets references unknown cell '" + p.node_name + "'");
            db_.addPin(db_.getCellId(p.node_name), net_id, p.offset[0], p.offset[1], p.direct ? std::string(1, p.direct) : std::string());
        }
    }
    void add_bookshelf_row(BookshelfParser::Row const& r) override {
        Row row; row.y = r.origin[1]; row.height = r.height; row.site_width = r.site_width; row.site_spacing = r.site_spacing; row.x_start = r.origin[0]; row.num_sites = r.site_num; row.x_end = row.x_start + row.num_sites * row.site_spacing; db_.addRow(row);
    }
    void set_bookshelf_node_position(std::string const& name, double x, double y, std::string const&, std::string const&, bool fixed) override { db_.setCellLocation(name, x, y, fixed); }
    void set_bookshelf_design(std::string&) override {}
    void bookshelf_end() override {}
private:
    PlacementDB& db_;
};
}
#endif

namespace {
std::string stripComment(const std::string& line) { auto p = line.find('#'); return p == std::string::npos ? line : line.substr(0, p); }
std::vector<std::string> split(const std::string& s) { std::istringstream iss(s); std::vector<std::string> out; for (std::string t; iss >> t;) out.push_back(t); return out; }
bool isNumberToken(const std::string& s) { if (s.empty()) return false; char* end = nullptr; std::strtod(s.c_str(), &end); return end == s.c_str() + s.size(); }
std::filesystem::path resolve(const std::filesystem::path& aux, const std::string& name) { auto p = std::filesystem::path(name); return p.is_absolute() ? p : aux.parent_path() / p; }

struct AuxFiles { std::string nodes, nets, pl, scl; };
AuxFiles parseAux(const std::string& aux_path) {
    std::ifstream in(aux_path);
    if (!in) throw std::runtime_error("Cannot open aux file: " + aux_path);
    std::string line; std::getline(in, line);
    for (char& ch : line) if (ch == ':' ) ch = ' ';
    auto toks = split(line);
    AuxFiles f;
    for (const auto& t : toks) {
        if (t.size() > 6 && t.substr(t.size()-6) == ".nodes") f.nodes = t;
        else if (t.size() > 5 && t.substr(t.size()-5) == ".nets") f.nets = t;
        else if (t.size() > 3 && t.substr(t.size()-3) == ".pl") f.pl = t;
        else if (t.size() > 4 && t.substr(t.size()-4) == ".scl") f.scl = t;
    }
    return f;
}

void parseNodes(PlacementDB& db, const std::filesystem::path& path) {
    std::ifstream in(path); if (!in) throw std::runtime_error("Cannot open nodes file: " + path.string());
    std::string line;
    while (std::getline(in, line)) {
        auto t = split(stripComment(line));
        if (t.size() >= 3 && t[0] != "UCLA" && isNumberToken(t[1]) && isNumberToken(t[2])) {
            bool terminal = t.size() >= 4 && t[3].find("terminal") != std::string::npos;
            db.addCell(t[0], std::stod(t[1]), std::stod(t[2]), terminal);
        }
    }
}
void parsePl(PlacementDB& db, const std::filesystem::path& path) {
    std::ifstream in(path); if (!in) throw std::runtime_error("Cannot open pl file: " + path.string());
    std::string line;
    while (std::getline(in, line)) {
        auto t = split(stripComment(line));
        if (t.size() >= 3 && t[0] != "UCLA") {
            bool fixed = line.find("/FIXED") != std::string::npos || line.find("/FIXED_NI") != std::string::npos;
            db.setCellLocation(t[0], std::stod(t[1]), std::stod(t[2]), fixed);
        }
    }
}
void parseNets(PlacementDB& db, const std::filesystem::path& path) {
    std::ifstream in(path); if (!in) throw std::runtime_error("Cannot open nets file: " + path.string());
    std::string line; int current_net = -1;
    while (std::getline(in, line)) {
        auto t = split(stripComment(line)); if (t.empty() || t[0] == "UCLA" || t[0].find(':') != std::string::npos) continue;
        if (t[0] == "NetDegree") { std::string name = t.size() >= 4 ? t[3] : "net_" + std::to_string(db.nets().size()); current_net = db.addNet(name); continue; }
        if (current_net >= 0 && t.size() >= 2) {
            if (!db.hasCell(t[0])) throw std::runtime_error(".nets references unknown cell '" + t[0] + "'");
            std::size_t offset_index = 2;
            if (offset_index < t.size() && t[offset_index] == ":") ++offset_index;
            double ox = offset_index < t.size() && isNumberToken(t[offset_index]) ? std::stod(t[offset_index]) : 0.0;
            double oy = offset_index + 1 < t.size() && isNumberToken(t[offset_index + 1]) ? std::stod(t[offset_index + 1]) : 0.0;
            db.addPin(db.getCellId(t[0]), current_net, ox, oy, t[1]);
        }
    }
}
void parseScl(PlacementDB& db, const std::filesystem::path& path) {
    std::ifstream in(path); if (!in) throw std::runtime_error("Cannot open scl file: " + path.string());
    std::string line; Row row; bool in_row = false;
    while (std::getline(in, line)) {
        auto t = split(stripComment(line)); if (t.empty()) continue;
        if (t[0] == "CoreRow") { row = Row{}; in_row = true; }
        else if (in_row && t[0] == "Coordinate" && t.size() >= 3) row.y = std::stod(t[2]);
        else if (in_row && t[0] == "Height" && t.size() >= 3) row.height = std::stod(t[2]);
        else if (in_row && t[0] == "Sitewidth" && t.size() >= 3) row.site_width = std::stod(t[2]);
        else if (in_row && t[0] == "Sitespacing" && t.size() >= 3) row.site_spacing = std::stod(t[2]);
        else if (in_row && t[0] == "SubrowOrigin" && t.size() >= 6) { row.x_start = std::stod(t[2]); row.num_sites = std::stoi(t[5]); }
        else if (in_row && t[0] == "End") { row.x_end = row.x_start + row.num_sites * row.site_spacing; db.addRow(row); in_row = false; }
    }
}
}

LimboBookshelfAdapter::LimboBookshelfAdapter(PlacementDB& db) : db_(db) {}

bool LimboBookshelfAdapter::read(const std::string& aux_path) {
#if defined(HAVE_LIMBO_BOOKSHELF)
    CallbackDB callback_db(db_);
    return BookshelfParser::read(callback_db, aux_path);
#else
    std::cerr << "Warning: built without Limbo headers; using a compatibility Bookshelf reader. Configure with -DLIMBO_ROOT=/path/to/limbo to call Limbo.\n";
    auto aux = std::filesystem::path(aux_path);
    auto files = parseAux(aux_path);
    parseNodes(db_, resolve(aux, files.nodes));
    parseNets(db_, resolve(aux, files.nets));
    parsePl(db_, resolve(aux, files.pl));
    if (!files.scl.empty()) parseScl(db_, resolve(aux, files.scl));
    return true;
#endif
}
