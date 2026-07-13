#pragma once

#include <string>
#include <vector>

/**
 * @brief Node type parsed from the Bookshelf .nodes file.
 */
enum class CellType {
    Standard,
    Terminal,
    TerminalNI
};


struct Cell {
    int id = -1;
    std::string name;
    double width = 0.0;
    double height = 0.0;
    double x = 0.0;
    double y = 0.0;

    CellType type = CellType::Standard;

// Fixed is a placement property, not a node-type property.
    bool is_fixed = false;

    std::vector<int> pin_ids;

    bool isTerminal() const {
        return type == CellType::Terminal ||
               type == CellType::TerminalNI;
    }

    bool isTerminalNI() const {
        return type == CellType::TerminalNI;
    }


    /** @brief Return true only for cells that may be moved by placement. */
    bool isMovable() const { return !is_fixed && type == CellType::Standard; }
};
