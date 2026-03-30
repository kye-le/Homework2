#include "geometry.h"

double cross(const Point& a, const Point& b) {
    return a.x * b.y - a.y * b.x;
}

double cross(const Point& a, const Point& b, const Point& c) {
    return cross(Point{b.x - a.x, b.y - a.y}, Point{c.x - a.x, c.y - a.y});
}

double dot(const Point& a, const Point& b) {
    return a.x * b.x + a.y * b.y;
}

double distance2(const Point& a, const Point& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

bool almost_equal(double a, double b, double eps) {
    return std::fabs(a - b) <= eps;
}

bool point_equal(const Point& a, const Point& b, double eps) {
    return almost_equal(a.x, b.x, eps) && almost_equal(a.y, b.y, eps);
}

int orientation(const Point& a, const Point& b, const Point& c) {
    const double val = cross(a, b, c);
    if (val > EPS) return 1;
    if (val < -EPS) return -1;
    return 0;
}

bool on_segment(const Point& a, const Point& b, const Point& p) {
    if (orientation(a, b, p) != 0) return false;
    return (std::min(a.x, b.x) - EPS <= p.x && p.x <= std::max(a.x, b.x) + EPS &&
            std::min(a.y, b.y) - EPS <= p.y && p.y <= std::max(a.y, b.y) + EPS);
}

bool segments_intersect(const Point& a, const Point& b, const Point& c, const Point& d) {
    const int o1 = orientation(a, b, c);
    const int o2 = orientation(a, b, d);
    const int o3 = orientation(c, d, a);
    const int o4 = orientation(c, d, b);

    if (o1 != o2 && o3 != o4) return true;

    if (o1 == 0 && on_segment(a, b, c)) return true;
    if (o2 == 0 && on_segment(a, b, d)) return true;
    if (o3 == 0 && on_segment(c, d, a)) return true;
    if (o4 == 0 && on_segment(c, d, b)) return true;

    return false;
}

double polygon_signed_area(const std::vector<Point>& pts) {
    const int n = static_cast<int>(pts.size());
    if (n < 3) return 0.0;

    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        const Point& p = pts[i];
        const Point& q = pts[(i + 1) % n];
        s += cross(p, q);
    }
    return 0.5 * s;
}

double triangle_signed_area2(const Point& a, const Point& b, const Point& c) {
    return cross(a, b) + cross(b, c) + cross(c, a);
}

Point area_preserving_point_baseline(const Point& A, const Point& B, const Point& C, const Point& D) {
    // Preserve local shoelace contribution:
    // cross(A,E) + cross(E,D) = S
    // where S = cross(A,B)+cross(B,C)+cross(C,D)
    //
    // Let E = M + t*n, where M is midpoint of AD and n is a unit normal to AD.
    // Then solve for t.
    //
    // This gives one valid area-preserving point; it is a baseline scaffold.
    const double S = cross(A, B) + cross(B, C) + cross(C, D);

    const Point M{(A.x + D.x) * 0.5, (A.y + D.y) * 0.5};
    const Point v{D.x - A.x, D.y - A.y};
    const double len = std::sqrt(v.x * v.x + v.y * v.y);

    if (len < EPS) {
        return M;
    }

    const Point n{-v.y / len, v.x / len};

    // S(E) = cross(A,E) + cross(E,D)
    // For E = M + t n:
    // S(E) = S(M) + t * (cross(A,n) + cross(n,D))
    const double SM = cross(A, M) + cross(M, D);
    const double coeff = cross(A, n) + cross(n, D);

    if (std::fabs(coeff) < EPS) {
        return M;
    }

    const double t = (S - SM) / coeff;
    return Point{M.x + t * n.x, M.y + t * n.y};
}

double local_displacement_proxy(const Point& A, const Point& B, const Point& C, const Point& D, const Point& E) {
    // Cheap baseline proxy:
    // sum absolute triangle areas between old chain and new chain.
    const double left = std::fabs(0.5 * triangle_signed_area2(A, B, E));
    const double mid  = std::fabs(0.5 * triangle_signed_area2(B, C, E));
    const double right= std::fabs(0.5 * triangle_signed_area2(C, D, E));
    return left + mid + right;
}