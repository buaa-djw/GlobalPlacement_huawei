#pragma once

#include <algorithm>

/**
 * @brief Axis-aligned rectangle used by placement geometry code.
 *
 * The box follows the half-open placement convention in spirit: positive area
 * overlap requires both width and height to be strictly positive, so touching
 * edges do not create bin membership.
 */
struct Box {
    double lx = 0.0;
    double ly = 0.0;
    double ux = 0.0;
    double uy = 0.0;

    /** @brief Return the x-span of the box. */
    double width() const { return ux - lx; }
    /** @brief Return the y-span of the box. */
    double height() const { return uy - ly; }
    /** @brief Return positive area, or zero for invalid/degenerate boxes. */
    double area() const { return valid() ? width() * height() : 0.0; }
    /** @brief A valid placement box must have strictly positive extents. */
    bool valid() const { return ux > lx && uy > ly; }
};

/**
 * @brief Compute strict positive-area overlap between two axis-aligned boxes.
 *
 * The overlap length in each dimension is the right/top minimum minus the
 * left/bottom maximum.  If either length is non-positive, the rectangles only
 * touch or are disjoint, so the function returns zero and callers must not add
 * spatial-index membership.
 */
inline double rectangleOverlapArea(const Box& lhs, const Box& rhs) {
    const double ox = std::min(lhs.ux, rhs.ux) - std::max(lhs.lx, rhs.lx);
    const double oy = std::min(lhs.uy, rhs.uy) - std::max(lhs.ly, rhs.ly);
    return (ox > 0.0 && oy > 0.0) ? ox * oy : 0.0;
}
