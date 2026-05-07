#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>


#define L 5.0
#define D 1.0
#define H L+D
#define W L
#define Lx W
#define Ly H
#define Re 100.0
#define Nu 1.0
#define poisson_tol 1e-5
#define poisson_max_iter 100000
#define V0 Re*Nu/D
double max_error;
int iter;

double **u;
double **v;
double **p;
double **u_star;
double **v_star;
double **rhs;

double dx, dy, dt;

int max_iter, Nx, Ny;

double **allocate_2d(int nx, int ny)
{
    int i;
    double **phi;

    phi = (double **)malloc((ny) * sizeof(double *));
    if (phi == NULL) {
        std::cout << "ERROR: malloc failed for row pointers\n";
        exit(1);
    }

    for (i = 0; i < ny; i++) {
        phi[i] = (double *)malloc((nx) * sizeof(double));
        if (phi[i] == NULL) {
            std::cout <<"ERROR: malloc failed for row " << i <<"\n";
            exit(1);
        }
    }

    return phi;
}

void free_2d(double **phi, int nx)
{
    int i;

    if (phi == NULL) return;

    for (i = 0; i < nx; i++) {
        free(phi[i]);
    }
    free(phi);
}

void free_all(){
    free_2d(p, Nx+2);
    free_2d(rhs, Nx+2);
    free_2d(u, Nx+1);
    free_2d(u_star, Nx+1);
    free_2d(v, Nx+2);
    free_2d(v_star, Nx+2);
}

void apply_grid(double dx, double dy) {
    Nx = (int)(Lx/dx);
    Ny = (int)((Ly)/dy);
    p = allocate_2d(Nx+2, Ny+2);
    rhs = allocate_2d(Nx+2, Ny+2);
    u = allocate_2d(Nx+1,Ny+2);
    u_star = allocate_2d(Nx+1, Ny+2);
    v = allocate_2d(Nx+2,Ny+1);
    v_star = allocate_2d(Nx+2, Ny+1);   
    for (int i = 0; i<Nx+2; i++) {
        for (int j = 0; j<Ny+2; j++) {
            p[j][i] = 0;
            rhs[j][i] = 0;
        }
    }
    for (int i = 0; i<Nx+1; i++) {
        for (int j = 0; j<Ny+2; j++) {
            u[j][i] = 0;
            u_star[j][i] = 0;
        }
    }
    for (int i = 0; i<Nx+2; i++) {
        for (int j = 0; j<Ny+1; j++) {
            v[j][i] = 0;
            v_star[j][i] = 0;
        }
    }
}

void apply_boundary_conditions() {
    for (int i = 0; i < Nx+1; i++) {
        u[0][i] = -u[1][i];
        u[Ny+1][i] = -u[Ny][i];
    }

    for (int j = 0; j < Ny+1; j++) {
        v[j][0] = -v[j][1];
        v[j][Nx+1] = -v[j][Nx];
    }

    for (int i = 1; i < Nx+1; i++) {
        v[0][i] = 0.0;
        v[Ny][i] = 0.0;
    }

    for (int j = 1; j < Ny+1; j++) {
        u[j][0] = 0.0;
        u[j][Nx] = 0.0;
    }

    // Inlet Boundary Conditions
    for (int i = 0; i < (int) (D/dx)+1; i++) {
        u[Ny+1][i] = -u[Ny][i];
    }

    for (int i = 1; i < (int) (D/dx)+1; i++) {
        v[Ny][i] = V0;
        p[Ny+1][i] = p[Ny][i];
    }

    // Outlet Boundary Conditions
    for (int i = (int) ((L-D)/dx); i < Nx+1; i++) {
        u[Ny+1][i] = u[Ny][i];
    }

    for (int i = (int) ((L-D)/dx) + 1; i< Nx+1; i++) {
        v[Ny][i] = v[Ny-1][i];
        p[Ny+1][i] = -p[Ny][i];
    }

    // Inner Side Walls
    for (int j = (int) (D/dy); j<Ny+1; j++) {
        v[j][(int) (D/dx)+1] = -v[j][(int) (D/dx)];
        v[j][(int) ((L-D)/dx)] = -v[j][(int) ((L-D)/dx)+1];
    }

    for (int j = (int) (D/dy)+1; j<Ny+1; j++) {
        u[j][(int) (D/dx)] = 0.0;
        u[j][(int) ((L-D)/dx)] = 0.0;
    }

    // Inner Bottom Wall
    for (int i = (int) (D/dx); i<(int) ((L-D)/dx); i++) {
        u[(int) (D/dy)+1][i] = -u[(int) (D/dy)][i];
    }

    for (int i = (int) (D/dx)+1; i<(int) ((L-D)/dx); i++) {
        v[(int) (D/dy)][i] = 0.0;
    }

    // Inner Volume
    for (int i = (int) (D/dx)+2; i<(int) ((L-D)/dx)-1; i++){
        for (int j = (int) (D/dy)+2; j<Ny+1; j++) {
            u[j][i] = 0;
            v[j][i] = 0;
            p[j][i] = 0;
        }
    }
}

void predictor() {
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx; i++) {
            if ((i>(D/dx))&(i<((L-D)/dx))&(j>(D/dy)+1)){
                u_star[j][i] = 0.0;
            } else {
                double ue = (u[j][i]+u[j][i+1])/2.0;
                double uw = (u[j][i]+u[j][i-1])/2.0;
                double un = (u[j][i]+u[j+1][i])/2.0;
                double us = (u[j][i]+u[j-1][i])/2.0;

                double vn = (v[j][i] + v[j][i+1])/2.0;
                double vs = (v[j-1][i] + v[j-1][i+1])/2.0;

                double Re_local = std::sqrt(u[j][i]*u[j][i] + v[j][i]*v[j][i])*dx/Nu;
                double conv;
                if (Re_local>2) {
                    conv = (ue*std::min(ue,0.0) + uw*std::min(-uw,0.0))*dy + (un*std::min(vn,0.0) + us*std::min(-vs,0.0))*dx;
                } else {
                    conv = (ue*ue - uw*uw)*dy + (un*vn - us*vs)*dx;
                }

                double diff = Nu*((u[j][i+1] - 2.0*u[j][i] + u[j][i-1])*dy/dx + (u[j+1][i] - 2.0*u[j][i] + u[j-1][i])*dx/dy);

                u_star[j][i] = u[j][i] + dt * (-conv + diff) / (dx*dy);
            }
        }
    }

    for (int j = 1; j < Ny; j++) {
        for (int i = 1; i < Nx+1; i++) {
            if ((i>(D/dx)+1)&(i<((L-D)/dx)-1)&(j>(D/dx))) {
                v_star[j][i] = 0.0;
            } else {
                double ve = (v[j][i]+v[j][i+1])/2.0;
                double vw = (v[j][i]+v[j][i-1])/2.0;
                double vn = (v[j][i]+v[j+1][i])/2.0;
                double vs = (v[j][i]+v[j-1][i])/2.0;

                double ue = (u[j][i]+u[j+1][i])/2.0;
                double uw = (u[j][i-1]+u[j+1][i-1])/2.0;

                double Re_local = std::sqrt(u[j][i]*u[j][i] + v[j][i]*v[j][i])*dx/Nu;
                double conv;
                if (Re_local>2) {
                    conv = (ve*std::min(ue,0.0) + vw*std::min(-uw,0.0))*dy + (vn*std::min(vn,0.0) + vs*std::min(-vs,0.0))*dx;
                } else {
                    conv = (ve*ue - vw*uw)*dy + (vn*vn - vs*vs)*dx;
                }

                double diff = Nu*((v[j][i+1] - 2.0*v[j][i] + v[j][i-1])*dy/dx + (v[j+1][i] - 2.0*v[j][i] + v[j-1][i])*dx/dy);

                v_star[j][i] = v[j][i] + dt * (-conv + diff) / (dx*dy);
            }
        }
    }
}

void find_divergence() {
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            rhs[j][i] = ((u_star[j][i] - u_star[j][i-1]) / dx +
                         (v_star[j][i] - v_star[j-1][i]) / dy) / dt;
        }
    }
}

void solve_poisson() {
    for (iter = 0; iter < poisson_max_iter; iter++) {
        max_error = 0.0;
        for (int j = 1; j < Ny+1; j++) {
            for (int i = 1; i < Nx+1; i++) {
                if ((i>(D/dx)+1)&(i<((L-D)/dx)-1)&(j>(D/dy)+1)) {
                    p[j][i] = 0.0;
                } else {
                    double p_old = p[j][i];
                    p[j][i] = (
                        (p[j][i+1] + p[j][i-1]) * dy * dy +
                        (p[j+1][i] + p[j-1][i]) * dx * dx -
                        rhs[j][i] * dx * dx * dy * dy
                    ) / (2.0 * (dx * dx + dy * dy));
                    double err = fabs(p[j][i] - p_old);
                    if (err > max_error) max_error = err;
                }
            }
        }
        /*for (int i = 0; i < Nx+2; i++) {
            p[0][i]     = p[1][i];       // bottom
            p[Ny+1][i]  = p[Ny][i];      // top
        }
        for (int j = 0; j < Ny+2; j++) {
            p[j][0]     = p[j][1];       // left
            p[j][Nx+1]  = p[j][Nx];      // right
        }*/
        // Fix pressure reference
        //p[Ny/2][Nx/2] = 0.0;

        if (max_error < poisson_tol) break;
    }
}

void corrector(){
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx; i++) {
            u[j][i] = u_star[j][i] - dt * (p[j][i+1] - p[j][i]) / dx;
        }
    }

    for (int j = 1; j < Ny; j++) {
        for (int i = 1; i < Nx+1; i++) {
            v[j][i] = v_star[j][i] - dt * (p[j+1][i] - p[j][i]) / dy;
        }
    }
}

void write_vtk(int step) {
    char filename[64];
    snprintf(filename, sizeof(filename), "output_%04d.vtk", step);
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        perror("File opening failed");
        return;
    }

    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "Lid-driven cavity flow\nASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", Nx, Ny);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING %f %f 1\n", dx, dy);
    fprintf(fp, "POINT_DATA %d\n", Nx * Ny);

    // Pressure field
    fprintf(fp, "SCALARS pressure float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            fprintf(fp, "%f\n", p[j][i]);
        }
    }

    // Velocity magnitude
    fprintf(fp, "SCALARS velocity_magnitude float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double uc = 0.5 * (u[j][i] + u[j][i-1]);
            double vc = 0.5 * (v[j][i] + v[j-1][i]);
            fprintf(fp, "%f\n", sqrt(uc * uc + vc * vc));
        }
    }

    // Vorticity (z-component of curl)
    fprintf(fp, "SCALARS vorticity float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double dudy = (u[j+1][i] - u[j-1][i]) / (2.0 * dy);
            double dvdx = (v[j][i+1] - v[j][i-1]) / (2.0 * dx);
            double omega = dvdx - dudy;
            fprintf(fp, "%f\n", omega);
        }
    }

    // u component
    fprintf(fp, "SCALARS u float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double uc = 0.5 * (u[j][i] + u[j][i-1]);
            fprintf(fp, "%f\n", uc);
        }
    }

    // v component
    fprintf(fp, "SCALARS v float 1\nLOOKUP_TABLE default\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double vc = 0.5 * (v[j][i] + v[j-1][i]);
            fprintf(fp, "%f\n", vc);
        }
    }

    // Velocity vector
    fprintf(fp, "VECTORS velocity float\n");
    for (int j = 1; j < Ny+1; j++) {
        for (int i = 1; i < Nx+1; i++) {
            double uc = 0.5 * (u[j][i] + u[j][i-1]);
            double vc = 0.5 * (v[j][i] + v[j-1][i]);
            fprintf(fp, "%f %f 0.0\n", uc, vc);
        }
    }

    fclose(fp);
}


int main(int argc, char *argv[]) {
    double t_max;
    if ((argc==4)&&(atof(argv[2])!=0.0)&&(atof(argv[3])!=0.0)){
        std::string dx_str = argv[1];
        if (dx_str.find("/") != std::string::npos) {
            int pos = dx_str.find("/");
            double numerator = std::stof(dx_str.substr(0,pos));
            double denominator = std::stof(dx_str.substr(pos+1));
            dx = numerator/denominator;
            dy = dx;
        } else {
            dx = atof(argv[1]);
            dy = dx;
        }
        dt = atof(argv[2]);
        t_max = atof(argv[3]);
    }
    else {
        dx = 1.0/16.0;
        dy = 1.0/16.0;
        dt = 1.0/1000.0;
        t_max = 10.0;
    }
    max_iter = (int) (t_max/dt);

    apply_grid(dx, dy);

    std::cout <<"dx = "<< dx <<",dy = " << dy <<", dt = "<< dt <<", max_iter = "<< max_iter<<"\n";

    for (int step = 0; step <= max_iter; step++) {
        apply_boundary_conditions();
        predictor();
        find_divergence();
        solve_poisson();
        corrector();

        printf("Step=%d iter=%d maxErr=%e u=%.4f v=%.4f p=%.4f\n",step, iter,
               max_error, u[(int)(D/(2*dy))][Nx/2], v[Ny/2][(int)(D/(2*dx))],p[(int)(D/(2*dy))][Nx/2]);
        if (step % 500 == 0) write_vtk(step);
    }

    free_all();
    std::cout << "Simulation complete.\n";
    return 0;
}
