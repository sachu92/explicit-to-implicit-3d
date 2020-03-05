/*
 *  Explicit to implicit reconstruction
 *  Date : 20 Feb 2020
 *  Author : Sachin Krishnan T V (sachu92@gmail.com)
 */

#ifndef __RECONSTRUCT_HPP__
#define __RECONSTRUCT_HPP__
#include <vector>

#define EPS 0.1

int mesh_size;
double xlo, xhi, ylo, yhi, zlo, zhi;
double xbinsize, ybinsize, zbinsize;

std::vector < std::vector < double > > pc_coord;
std::vector < std::vector < double > > pc_norm;

std::vector < double > rbf_weight;

std::vector < double > mesh_ls;

void readPointCloud(char *);
double distance(double, double, double, double, double, double);
double triharmonic_kernel(double);
void evaluateSDF();
void reconstructMesh();
void outputMesh();

#endif
