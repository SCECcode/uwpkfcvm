#include <stdio.h>
#include <math.h>
#include <stdbool.h>

/* ------------------ basic vector ------------------ */
typedef struct {
    double x, y, z;
} Vec3;

Vec3 vsub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 vcross(Vec3 a, Vec3 b) {
    return (Vec3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

double vdot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

/* ------------------ barycentric ------------------ */
bool barycentric_tet(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3,
                     Vec3 p,
                     double *l0, double *l1, double *l2, double *l3)
{
    Vec3 a = vsub(p1, p0);
    Vec3 b = vsub(p2, p0);
    Vec3 c = vsub(p3, p0);
    Vec3 ap = vsub(p, p0);

    double denom = vdot(a, vcross(b, c));
    if (fabs(denom) < 1e-20) return false;

    double w1 = vdot(ap, vcross(b, c)) / denom;
    double w2 = vdot(a, vcross(ap, c)) / denom;
    double w3 = vdot(a, vcross(b, ap)) / denom;
    double w0 = 1.0 - w1 - w2 - w3;

    *l0 = w0; *l1 = w1; *l2 = w2; *l3 = w3;
    return true;
}

bool point_in_tet(Vec3 p0, Vec3 p1, Vec3 p2, Vec3 p3,
                  Vec3 p,
                  double *l0, double *l1, double *l2, double *l3)
{
    if (!barycentric_tet(p0,p1,p2,p3,p,l0,l1,l2,l3)) return false;

    return (*l0 >= -1e-12 && *l1 >= -1e-12 &&
            *l2 >= -1e-12 && *l3 >= -1e-12);
}

/* ------------------ tetra interpolation ------------------ */
double interp_tets(Vec3 q, Vec3 p[8], double val[8])
{
    int tets[6][4] = {
        {0,1,3,7},
        {0,3,2,7},
        {0,2,6,7},
        {0,6,4,7},
        {0,4,5,7},
        {0,5,1,7}
    };

    for (int t = 0; t < 6; t++) {
        int i0 = tets[t][0];
        int i1 = tets[t][1];
        int i2 = tets[t][2];
        int i3 = tets[t][3];

        double l0,l1,l2,l3;
        if (point_in_tet(p[i0],p[i1],p[i2],p[i3],
                         q,&l0,&l1,&l2,&l3)) {

            return l0*val[i0] +
                   l1*val[i1] +
                   l2*val[i2] +
                   l3*val[i3];
        }
    }

    return NAN;
}

/* ------------------ main test ------------------ */
int main()
{
    /* define cube corners */
    Vec3 p[8] = {
        {0,0,0}, // 0
        {1,0,0}, // 1
        {0,1,0}, // 2
        {1,1,0}, // 3
        {0,0,1}, // 4
        {1,0,1}, // 5
        {0,1,1}, // 6
        {1,1,1}  // 7
    };

    /* define scalar field (example: f(x,y,z) = x + y + z) */
    double val[8];
    for (int i = 0; i < 8; i++) {
        val[i] = p[i].x + p[i].y + p[i].z;
    }

    /* query point inside cube */
    Vec3 q = {0.3, 0.2, 0.4};

    double result = interp_tets(q, p, val);

    printf("Interpolated value: %f\n", result);
    printf("Expected value:     %f\n", q.x + q.y + q.z);

    return 0;
}

