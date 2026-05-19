#include <math.h>
#include <stdbool.h>

typedef struct {
    double x, y, z;
} KDV3_simple;

static inline KDV3_simple vsub(KDV3_simple a, KDV3_simple b) {
    KDV3_simple r = {a.x - b.x, a.y - b.y, a.z - b.z};
    return r;
}

static inline KDV3_simple vcross(KDV3_simple a, KDV3_simple b) {
    KDV3_simple r = {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
    return r;
}

static inline double vdot(KDV3_simple a, KDV3_simple b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static bool barycentric_tet(KDV3_simple p0, KDV3_simple p1, KDV3_simple p2, KDV3_simple p3,
                            KDV3_simple p, double *l0, double *l1, double *l2, double *l3)
{
    KDV3_simple a = vsub(p1, p0);
    KDV3_simple b = vsub(p2, p0);
    KDV3_simple c = vsub(p3, p0);
    KDV3_simple ap = vsub(p,  p0);

    double denom = vdot(a, vcross(b, c));
    if (fabs(denom) < 1e-20) {
        return false; // degenerate tetra
    }

    double w1 = vdot(ap, vcross(b, c)) / denom;
    double w2 = vdot(a,  vcross(ap, c)) / denom;
    double w3 = vdot(a,  vcross(b, ap)) / denom;
    double w0 = 1.0 - w1 - w2 - w3;

    if (l0) *l0 = w0;
    if (l1) *l1 = w1;
    if (l2) *l2 = w2;
    if (l3) *l3 = w3;

    return true;
}

static bool point_in_tet(KDV3_simple p0, KDV3_simple p1, KDV3_simple p2, KDV3_simple p3,
                         KDV3_simple p, double eps,
                         double *l0, double *l1, double *l2, double *l3)
{
    if (!barycentric_tet(p0, p1, p2, p3, p, l0, l1, l2, l3)) {
        return false;
    }

    return (*l0 >= -eps && *l1 >= -eps && *l2 >= -eps && *l3 >= -eps);
}

static double tet_interp_scalar(KDV3_simple q,
                                KDV3_simple p0, KDV3_simple p1, KDV3_simple p2, KDV3_simple p3,
                                double v0, double v1, double v2, double v3,
                                bool *ok)
{
    double l0, l1, l2, l3;
    if (!point_in_tet(p0, p1, p2, p3, q, 1e-12, &l0, &l1, &l2, &l3)) {
        if (ok) *ok = false;
        return NAN;
    }

    if (ok) *ok = true;
    return l0 * v0 + l1 * v1 + l2 * v2 + l3 * v3;
}

/* Try the 6 tetrahedra of a hexahedral cell */
static double interp_cell_tets(KDV3_simple q,
                               KDV3_simple p[8],
                               double v[8],
                               bool *ok)
{
    static const int tets[6][4] = {
        {0, 1, 3, 7},
        {0, 3, 2, 7},
        {0, 2, 6, 7},
        {0, 6, 4, 7},
        {0, 4, 5, 7},
        {0, 5, 1, 7}
    };

    for (int t = 0; t < 6; ++t) {
        int i0 = tets[t][0];
        int i1 = tets[t][1];
        int i2 = tets[t][2];
        int i3 = tets[t][3];

        double out = tet_interp_scalar(q,
                                       p[i0], p[i1], p[i2], p[i3],
                                       v[i0], v[i1], v[i2], v[i3],
                                       ok);
        if (ok && *ok) {
            return out;
        }
    }

    if (ok) *ok = false;
    return NAN;
}


void uwpkfcvm_read_interp_properties(uwpkfcvm_model_t *model, int lldindex,
                                     uwpkfcvm_properties_t *data,
                                     double lat, double lon, double depth)
{
    int nx = model->nx;
    int ny = model->ny;
    int nz = model->nz;   /* not model->ny */
    int sz = model->pnts_size;

    KDVec3 *xyz = model->v3pnts;
    KDVec3 query_xyz;

    lld_to_xyz(&query_xyz, lat, lon, depth, -1);

    int xidx, yidx, zidx;
    lldindex_to_idx(lldindex, nx, ny, &xidx, &yidx, &zidx);

    int offset[8];
    offset[0] = idx_to_lldindex(nx, ny, xidx,     yidx,     zidx);
    offset[1] = idx_to_lldindex(nx, ny, xidx + 1, yidx,     zidx);
    offset[2] = idx_to_lldindex(nx, ny, xidx,     yidx + 1, zidx);
    offset[3] = idx_to_lldindex(nx, ny, xidx + 1, yidx + 1, zidx);
    offset[4] = idx_to_lldindex(nx, ny, xidx,     yidx,     zidx + 1);
    offset[5] = idx_to_lldindex(nx, ny, xidx + 1, yidx,     zidx + 1);
    offset[6] = idx_to_lldindex(nx, ny, xidx,     yidx + 1, zidx + 1);
    offset[7] = idx_to_lldindex(nx, ny, xidx + 1, yidx + 1, zidx + 1);

    KDVec3 p[8];
    double vs[8];
    double vp[8];

    for (int i = 0; i < 8; ++i) {
        p[i]  = *find_xyz_by_lldindex(xyz, sz, offset[i]);
        vs[i] = vs_by_offset(model, offset[i]);
        vp[i] = vp_by_offset(model, offset[i]);
    }

    bool ok_vs = false;
    bool ok_vp = false;

    data->vs = interp_cell_tets(query_xyz, p, vs, &ok_vs);
    data->vp = interp_cell_tets(query_xyz, p, vp, &ok_vp);

    if (!ok_vs || !ok_vp) {
        /* fallback if query is outside the cell or the cell is degenerate */
        data->vs = NAN;
        data->vp = NAN;
    }

    if (data->vp > 0.0) {
        data->rho = uwpkfcvm_calculate_density(data->vp);
    }
}
