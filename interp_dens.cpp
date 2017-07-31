
/*
#
# interp_dens.cpp
# by Mauro Alberti - www.malg.eu, alberti.m65@gmail.com
# 2011-03-13->05-04
# vers. 0.8
#
# This program or module is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 2 of the License, or
# version 3 of the License, or (at your option) any later version. It is
# provided for educational purposes and is distributed in the hope that
# it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
# the GNU General Public License for more details.
#
*/

#include <stdio.h>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>



using namespace std;


class Point3d
{
    public:
        //data members
        double x, y, z;
        // Constructor uses default arguments to allow calling with zero, one,
        // or two values.
        Point3d(double xval = 0.0, double yval = 0.0, double zval = 0.0)
        {
                x = xval;
                y = yval;
                z = zval;
        }
        // Distance to another point.
        double distance(Point3d other)
        {
                double xd = x - other.x;
                double yd = y - other.y;
                double zd = z - other.z;
                return sqrt(xd*xd + yd*yd + zd*zd);
        }
};



int isfloat (string instring)
{
    for (int i = 0; i < instring.size(); i++)
    {
        if (!((instring[i]>= '0' && instring[i] <='9') or (instring[i] == '.')))
        {
            return 0;
        }
    }
    return 1;
}


double calc_normkern_val(double pts_distance, double kernel_bandwidth)
{
    double t = pts_distance / kernel_bandwidth;
    return exp(-t*t/2);
}

double calc_quartkern_val(double pts_distance, double kernel_bandwidth)
{
    double t = pts_distance / kernel_bandwidth;
    t = 1.0 - (t*t);
    return t*t;
}

int main ()
{

    double const Pi=4*atan(1);

    char comma;

    cout << "\n\nDensityInterpolation\n";
    cout << "by M. Alberti - www.malg.eu\n";
    cout << "v. 1.0 - 2011-04\n\n\n";


    // input of parameter file name and file opening
    cout << "Enter name of input parameter file: ";
  	string param_filename;
    cin >> param_filename;

    ifstream param_file;
    param_file.open(param_filename.c_str(), ios::binary);

    if (param_file.fail())
    {
    cout << "\n\nUnable to open parameter file '" << param_filename << "'\n";
    return 1;
    }

    // start time
    clock_t start_time = clock();


    // parameter reading

    string data_filename_in, grid_density_flnm_in, cell_size_str_in, kernel_functions_in, kernel_bandwidth_in;
    string data_filename, grid_density_flnm, cell_size_str, kernel_functions_str, kernel_bandwidth_str;

    cout << "\nParameters\n";
    getline(param_file, data_filename_in); istringstream df(data_filename_in); df >> data_filename; cout << " - 3D point data (input): " << data_filename << "\n";
    getline(param_file, grid_density_flnm_in); istringstream gd(grid_density_flnm_in); gd >> grid_density_flnm; cout << " - grid of point density (output): " << grid_density_flnm << "\n";
    getline(param_file, cell_size_str_in); istringstream cs(cell_size_str_in); cs >> cell_size_str; cout << " - cell size: " << cell_size_str << "\n";
    getline(param_file, kernel_functions_in); istringstream kf(kernel_functions_in); kf >> kernel_functions_str; cout << " - kernel functions: " << kernel_functions_str << "\n";
    getline(param_file, kernel_bandwidth_in); istringstream kb(kernel_bandwidth_in); kb >> kernel_bandwidth_str; cout << " - kernel bandwidth: " << kernel_bandwidth_str << "\n";


    // test cell size value being number
    if (!isfloat(cell_size_str))
    {
        cout << "\n\nInput cell size is not number\n";
        return 1;
    }
    double cell_size;
	istringstream instr(cell_size_str);
	instr >> cell_size;

    // test and process chosen kernel functions
    string kernel_functions;
	istringstream instr_kf(kernel_functions_str);
	instr_kf >> kernel_functions;
	bool normal_kernel = false, quartic_kernel = false;
    for (int i = 0; i < kernel_functions.size(); i++)
    {
        if (kernel_functions[i] == 'n')
        {
            normal_kernel = true;
        }
        else if (kernel_functions[i] == 'q')
        {
            quartic_kernel = true;
        }
        else
        {
            cout << "\n\nKernel functions are not correctly defined (n and/or q)\n";
            return 1;
        }
    }


    // test kernel bandwidth value being number
    if (!isfloat(kernel_bandwidth_str))
    {
        cout << "\n\nInput kernel bandwidth is not number\n";
        return 1;
    }
    double kernel_bandwidth;
	istringstream instr_kb(kernel_bandwidth_str);
	instr_kb >> kernel_bandwidth;


    // input file processing
    ifstream infile;
    infile.open(data_filename.c_str(), ios::binary);

    if (infile.fail())
    {
    cout << "\n\nUnable to open input file '" << data_filename << "'\n";
    return 1;
    }


    // output file processing
    FILE * outgridfile;
    outgridfile = fopen (grid_density_flnm.c_str(),"w");


    //
    // calculations begin
    //

    string rec_line;

    //  list of raw data strings
    list<string> rawdata_list;

    // read file header
    getline(infile, rec_line);

    while (!infile.eof())
    {
        getline(infile, rec_line);
        if (rec_line.size() > 0)
        {
            rawdata_list.push_back(rec_line);
        }
    }

    infile.close();

    int num_recs = rawdata_list.size();

    list<string>::iterator string_pos; // string list iterator

    vector<double> x(num_recs), y(num_recs), z(num_recs);

    int rec_ndx = 0;

    for (string_pos=rawdata_list.begin(); string_pos!=rawdata_list.end(); string_pos++)
    {

        istringstream instr(*string_pos);
        instr >> x[rec_ndx] >> comma >> y[rec_ndx] >> comma >> z[rec_ndx];
        rec_ndx++;

    }

    // spatial range: min and max of x, y and z
    double x_min = *min_element( x.begin(), x.end() );
    double x_max = *max_element( x.begin(), x.end() );
    double y_min = *min_element( y.begin(), y.end() );
    double y_max = *max_element( y.begin(), y.end() );
    double z_min = *min_element( z.begin(), z.end() );
    double z_max = *max_element( z.begin(), z.end() );

    cout << "Spatial range of data:\n";
    cout << "  x_min: " << x_min << " x_max: " << x_max << "\n";
    cout << "  y_min: " << y_min << " y_max: " << y_max << "\n";
    cout << "  z_min: " << z_min << " z_max: " << z_max << "\n\n";

    float x_range = x_max -  x_min;
    float y_range = y_max -  y_min;
    float z_range = z_max -  z_min;

    int num_x_lines = int(x_range/cell_size) + 1;
    int num_y_lines = int(y_range/cell_size) + 1;
    int num_z_lines = int(z_range/cell_size) + 1;

    cout << "Columns (x) = " << num_x_lines << "; rows (y) = " << num_y_lines << "; z lines = " << num_z_lines << "\n";


    // definies 3d array storing values for normal distribution kernel
    // modified from http://www.cplusplus.com/forum/articles/7459/ (consulted 2011-05-03)
    vector<vector<vector<double> > > normal_kernel_values;
    // Set up sizes. (num_y_lines x num_x_lines)
    normal_kernel_values.resize(num_y_lines);
    for (int i = 0; i < num_y_lines; i++)
    {
        normal_kernel_values[i].resize(num_x_lines);
        for (int j = 0; j < num_x_lines; j++)
          normal_kernel_values[i][j].resize(num_z_lines);
    }

    // defines 3d array storing values for quartic distribution kernel
    // modified from http://www.cplusplus.com/forum/articles/7459/ (consulted 2011-05-03)
    vector<vector<vector<double> > > quartic_kernel_values;
    // Set up sizes. (num_y_lines x num_x_lines)
    quartic_kernel_values.resize(num_y_lines);
    for (int i = 0; i < num_y_lines; i++)
    {
        quartic_kernel_values[i].resize(num_x_lines);
        for (int j = 0; j < num_x_lines; j++)
          quartic_kernel_values[i][j].resize(num_z_lines);
    }

    // calculation of grid indices for 3D points - i downwards, j to the right, k towards the observer
   vector<double> rec_ndx_i(num_recs), rec_ndx_j(num_recs), rec_ndx_k(num_recs);

    for (int ndx = 0; ndx < num_recs; ndx++)
    {
        rec_ndx_i[ndx] = (y_max - y[ndx])/cell_size;
        rec_ndx_j[ndx] = (x[ndx]- x_min)/cell_size;
        rec_ndx_k[ndx] = (z[ndx]- z_min)/cell_size;
    }


    cout << "\nCalculating .... Please wait\n";

    fprintf (outgridfile, "# vtk DataFile Version 3.0\n");
    fprintf (outgridfile, "Density interpolation\n\n");
    fprintf (outgridfile, "ASCII\n");
    fprintf (outgridfile, "DATASET STRUCTURED_POINTS\n");
    fprintf (outgridfile, "DIMENSIONS %i %i %i\n", num_x_lines, num_y_lines, num_z_lines );
    fprintf (outgridfile, "ORIGIN %f %f %f\n", x_min+(cell_size/2.0), y_min+(cell_size/2.0), z_min+(cell_size/2.0) );
    fprintf (outgridfile, "SPACING %f %f %f\n", cell_size, cell_size, cell_size );
    fprintf (outgridfile, "POINT_DATA %i\n", num_x_lines*num_y_lines*num_z_lines );


    for (int i = 0; i < num_y_lines; i++)
    {
        for (int j = 0; j < num_x_lines; j++)
        {
            for (int k = 0; k < num_z_lines; k++)
            {

                Point3d grid_pt = Point3d(i+0.5, j+0.5, k+0.5);
                double normkern_cellvalue = 0.0, quartkern_cellvalue = 0.0;
                for (int ndx = 0; ndx < num_recs; ndx++)
                {
                    Point3d curr_pt = Point3d(rec_ndx_i[ndx], rec_ndx_j[ndx], rec_ndx_k[ndx]);
                    double pts_distance = grid_pt.distance(curr_pt);

                    if (normal_kernel)
                    {
                        normkern_cellvalue += calc_normkern_val(pts_distance, kernel_bandwidth);
                    }
                    if (quartic_kernel)
                    {
                        if (pts_distance > kernel_bandwidth) continue;
                        quartkern_cellvalue += calc_quartkern_val(pts_distance, kernel_bandwidth);

                    }


                }
                if (normal_kernel) normal_kernel_values[i][j][k] = normkern_cellvalue/(2*Pi*kernel_bandwidth*kernel_bandwidth);
                if (quartic_kernel) quartic_kernel_values[i][j][k] = quartkern_cellvalue*3/(Pi*kernel_bandwidth*kernel_bandwidth);

                //
            }
        }
    }


    if (normal_kernel)
    {
        fprintf (outgridfile, "SCALARS density_normal float 1\n");
        fprintf (outgridfile, "LOOKUP_TABLE default\n");
        for (int k = 0; k < num_z_lines; k++)
        {
            for (int i = num_y_lines - 1; i >= 0; i--)
            {
                for (int j = 0; j < num_x_lines; j++)
                {
                    fprintf (outgridfile, "%f\n", normal_kernel_values[i][j][k]);
                }
            }
        }
    }

    if (quartic_kernel)
    {
        fprintf (outgridfile, "\n\nSCALARS density_quartic float 1\n");
        fprintf (outgridfile, "LOOKUP_TABLE default\n");
        for (int k = 0; k < num_z_lines; k++)
        {
            for (int i = num_y_lines - 1; i >= 0; i--)
            {
                for (int j = 0; j < num_x_lines; j++)
                {
                    fprintf (outgridfile, "%f\n", quartic_kernel_values[i][j][k]);
                }
            }
        }
    }

    fclose (outgridfile);

    // end time
    clock_t end_time = clock();
    float diff_time = ((float)end_time - (float)start_time)/1000.0;  // run time

    printf ("\n\nProcessing completed in %.2lf seconds\n\n", diff_time );

    return 0;
}




