#ifndef TYPEEXG_DLIB_ARMA_H
#define TYPEEXG_DLIB_ARMA_H

// Copyright (C) 2017 Kyaw Kyaw Htike @ Ali Abdul Ghafur. All rights reserved.

#include <vector>
#include "armadillo"
#include "dlib/image_transforms.h"
#include "dlib/matrix.h"
#include "dlib/array2d.h"
#include "dlib/array.h"

// RGB color image
// This function makes a deep copy.
void dlib2arma(const dlib::array2d<dlib::rgb_pixel> &mat_in, arma::Cube<unsigned char> &mat_out)
{

	int nrows = dlib::num_rows(mat_in);
	int ncols = dlib::num_columns(mat_in);
	const int nchannels = 3;

	mat_out.set_size(nrows, ncols, nchannels);

	unsigned char * ptr = mat_out.memptr();
	unsigned long count = 0;
	
	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].red;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].green;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].blue;
}

// BGR color image
void dlib2arma(const dlib::array2d<dlib::bgr_pixel> &mat_in, arma::Cube<unsigned char> &mat_out)
{

	int nrows = dlib::num_rows(mat_in);
	int ncols = dlib::num_columns(mat_in);
	const int nchannels = 3;

	mat_out.set_size(nrows, ncols, nchannels);

	unsigned char * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].blue;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].green;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].red;
}

// HSI color image
void dlib2arma(const dlib::array2d<dlib::hsi_pixel> &mat_in, arma::Cube<unsigned char> &mat_out)
{

	int nrows = dlib::num_rows(mat_in);
	int ncols = dlib::num_columns(mat_in);
	const int nchannels = 3;

	mat_out.set_size(nrows, ncols, nchannels);

	unsigned char * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].h;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].s;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].i;
}

// LAB color image
void dlib2arma(const dlib::array2d<dlib::lab_pixel> &mat_in, arma::Cube<unsigned char> &mat_out)
{

	int nrows = dlib::num_rows(mat_in);
	int ncols = dlib::num_columns(mat_in);
	const int nchannels = 3;

	mat_out.set_size(nrows, ncols, nchannels);

	unsigned char * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].l;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].a;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].b;
}

// RGBA color image
void dlib2arma(const dlib::array2d<dlib::rgb_alpha_pixel> &mat_in, arma::Cube<unsigned char> &mat_out)
{

	int nrows = dlib::num_rows(mat_in);
	int ncols = dlib::num_columns(mat_in);
	const int nchannels = 4;

	mat_out.set_size(nrows, ncols, nchannels);

	unsigned char * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].red;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].green;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].blue;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j].alpha;
}

// T can be unsigned char, unsigned short, unsigned int, unsigned long,
// char, signed char, short, int, long
template <typename T>
void dlib2arma(const dlib::array2d<T> &mat_in, arma::Mat<T> &mat_out)
{
	int nrows = dlib::num_rows(mat_in);
	int ncols = dlib::num_columns(mat_in);

	mat_out.set_size(nrows, ncols);

	T * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in[i][j];
}

// T should be generic C++ native types
// the 2nd arg mat_out will always be a 2D matrix with single channel
// this function makes a deep copy
template <typename T>
void dlib2arma(const dlib::matrix<T> &mat_in, arma::Mat<T> &mat_out)
{
	int nrows = mat_in.nr();
	int ncols = mat_in.nc();

	mat_out.set_size(nrows, ncols);

	T * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			ptr[count++] = mat_in(i,j);
}

// T should be generic C++ native types
// the 2nd arg mat_out will always be a column vector
// this function makes a deep copy
template <typename T>
void dlib2arma(const dlib::matrix<T,0,1> &mat_in, arma::Col<T> &mat_out)
{
	int nrows = mat_in.nr();
	//int ncols = mat_in.nc(); // should be always 1

	mat_out.set_size(nrows);

	T * ptr = mat_out.memptr();

	for (unsigned int i = 0; i < nrows; i++)
			ptr[i] = mat_in(i);
}

// T should be generic C++ native types
// the 2nd arg mat_out will always be a row vector
// this function makes a deep copy
template <typename T>
void dlib2arma(const dlib::matrix<T, 1, 0> &mat_in, arma::Row<T> &mat_out)
{
	//int nrows = mat_in.nr(); // should be always 1
	int ncols = mat_in.nc();

	mat_out.set_size(ncols);

	T * ptr = mat_out.memptr();

	for (unsigned int i = 0; i < ncols; i++)
		ptr[i] = mat_in(i);
}

// T should be generic C++ native types
// mat_in (the 1st arg) should be either only one row or one column (i.e. similar to a vector)
// this function makes a deep copy
template <typename T, int nchannels>
void dlib2arma(const dlib::array2d<dlib::matrix<T, nchannels, 1L>> &mat_in, arma::Cube<T> &mat_out)
{
	int nr = mat_in.nr();
	int nc = mat_in.nc();
	//int nchannels = mat_in[0][0].nr() * mat_in[0][0].nc();

	mat_out.set_size(nr, nc, nchannels);

	T * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int k = 0; k < nchannels; k++)
		for (int j = 0; j < nc; j++)
			for (int i = 0; i < nr; i++)
				ptr[count++] = mat_in[i][j](k);
}

// T should be generic C++ native types
// this function makes a deep copy
template <typename T, int nchannels>
void dlib2arma(const dlib::array<dlib::array2d<T>> &mat_in, arma::Cube<T> &mat_out)
{
	int nr = mat_in[0].nr();
	int nc = mat_in[0].nc();
	//int nchannels = mat_in.size();

	mat_out.set_size(nr, nc, nchannels);

	T * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int k = 0; k < nchannels; k++)
		for (int j = 0; j < nc; j++)
			for (int i = 0; i < nr; i++)
				ptr[count++] = mat_in[k][i][j];

}

// T should be generic C++ native types
template <typename T, int nchannels>
void dlib2arma(const dlib::array<dlib::matrix<T>> &mat_in, arma::Cube<T> &mat_out)
{
	int nr = mat_in[0].nr();
	int nc = mat_in[0].nc();
	//int nchannels = mat_in.size();

	mat_out.set_size(nr, nc, nchannels);

	dlib::matrix<T> temp;

	T * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int k = 0; k < nchannels; k++)
	{
		temp = mat_in[k];
		for (int j = 0; j < nc; j++)
			for (int i = 0; i < nr; i++)
				ptr[count++] = temp(i,j);
	}
}

// T should be generic C++ native types
template <typename T>
void dlib2arma(const std::vector<dlib::matrix<T>> &mat_in, arma::Cube<T> &mat_out)
{
	int nr = mat_in[0].nr();
	int nc = mat_in[0].nc();
	int nchannels = mat_in.size();

	mat_out.set_size(nr, nc, nchannels);

	dlib::matrix<T> t;

	T * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int k = 0; k < nchannels; k++)
	{
		t = mat_in[k];
		for (int j = 0; j < nc; j++)
			for (int i = 0; i < nr; i++)
				ptr[count++] = t(i, j);
	}
}

// T should be generic C++ native types
template <typename T, int nrows, int ncols>
void dlib2arma(const std::vector<dlib::matrix<T, nrows, ncols>> &mat_in, arma::Cube<T> &mat_out)
{
	//int nr = mat_in[0].nr();
	//int nc = mat_in[0].nc();
	int nr = nrows; 
	int nc = ncols;
	int nchannels = mat_in.size();

	mat_out.set_size(nr, nc, nchannels);

	dlib::matrix<T, nrows, ncols> temp;

	T * ptr = mat_out.memptr();
	unsigned long count = 0;

	for (int k = 0; k < nchannels; k++)
	{
		temp = mat_in[k];
		for (int j = 0; j < nc; j++)
			for (int i = 0; i < nr; i++)
				ptr[count++] = temp(i, j);
	}
}


// T should be generic C++ native types
// this is obviously only for 2D matrices (1 channel)
// makes deep copy
template <typename T>
void arma2dlib(arma::Mat<T> mat_in, dlib::matrix<T> &mat_out)
{
	mat_out = dlib::mat(mat_in);
}


// T can be unsigned char, unsigned short, unsigned int, unsigned long,
// char, signed char, short, int, long
// output mat_in can be used as grayscale image or for any other purpose
// makes deep copy
template <typename T>
void arma2dlib(const arma::Mat<T> &mat_in, dlib::array2d<T> &mat_out)
{
	int nrows = mat_in.n_rows;
	int ncols = mat_in.n_cols;

	mat_out.set_size(nrows, ncols);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j] = mat_in.at(i,j);		
}

// T should be generic C++ native types
// the 2nd arg mat_out will always be a column vector
// this function makes a deep copy
template <typename T>
void arma2dlib(const arma::Col<T> &mat_in, dlib::matrix<T, 0, 1> &mat_out)
{
	
	int nrows = mat_in.n_rows;
	//int ncols = mat_in.n_cols; // should be always 1

	mat_out.set_size(nrows);

	for (unsigned int i = 0; i < nrows; i++)
		mat_out(i) = mat_in.at(i);
}

// T should be generic C++ native types
// the 2nd arg mat_out will always be a row vector
// this function makes a deep copy
template <typename T>
void arma2dlib(const arma::Row<T> &mat_in, dlib::matrix<T, 1, 0> &mat_out)
{
	//int nrows = mat_in.n_rows; // should be always 1
	int ncols = mat_in.n_cols;

	mat_out.set_size(ncols);

	for (unsigned int i = 0; i < ncols; i++)
		mat_out(i) = mat_in.at(i);
}

// T should be generic C++ native types
// mat_in (the 1st arg) should be either only one row or one column (i.e. similar to a vector)
// this function makes a deep copy
template <typename T, int nchannels>
void arma2dlib(const arma::Cube<T> &mat_in, dlib::array2d<dlib::matrix<T, nchannels, 1L>> &mat_out)
{
	int nr = mat_in.n_rows;
	int nc = mat_in.n_cols;
	//int nchannels = mat_in.n_slices;

	mat_out.set_size(nr, nc);

	for (int k = 0; k < nchannels; k++)
		for (int j = 0; j < nc; j++)
			for (int i = 0; i < nr; i++)
				mat_out[i][j](k) = mat_in.at(i,j,k);
}

// T should be generic C++ native types
// this function makes a deep copy
template <typename T, int nchannels>
void arma2dlib(const arma::Cube<T> &mat_in, dlib::array<dlib::array2d<T>> &mat_out)
{
	int nr = mat_in.n_rows;
	int nc = mat_in.n_cols;
	//int nchannels = mat_in.n_slices;

	mat_out.set_max_size(nchannels);
	mat_out.set_size(nchannels);

	for (int k = 0; k < nchannels; k++)
		mat_out[k].set_size(nr, nc);

	for (int k = 0; k < nchannels; k++)
		for (int j = 0; j < nc; j++)
			for (int i = 0; i < nr; i++)
				mat_out[k][i][j] = mat_in.at(i,j,k);
}

// T should be generic C++ native types
template <typename T, int nchannels>
void arma2dlib(const arma::Cube<T> &mat_in, dlib::array<dlib::matrix<T>> &mat_out)
{
	int nr = mat_in.n_rows;
	int nc = mat_in.n_cols;
	//int nchannels = mat_in.n_slices;

	mat_out.set_max_size(nchannels);
	mat_out.set_size(nchannels);

	dlib::matrix<T> temp(nr, nc);
	
	for (int k = 0; k < nchannels; k++)
	{
		for (int j = 0; j < nc; j++)
		{
			for (int i = 0; i < nr; i++)
				temp(i,j) = mat_in.at(i,j,k);
		}
		mat_out[k] = temp;
	}
}

// T should be generic C++ native types
template <typename T>
void arma2dlib(const arma::Cube<T> &mat_in, std::vector<dlib::matrix<T>> &mat_out)
{
	int nr = mat_in.n_rows;
	int nc = mat_in.n_cols;
	int nchannels = mat_in.n_slices;

	mat_out.resize(nchannels);

	dlib::matrix<T> temp(nr, nc);

	for (int k = 0; k < nchannels; k++)
	{
		for (int j = 0; j < nc; j++)
		{
			for (int i = 0; i < nr; i++)
				temp(i, j) = mat_in.at(i,j,k);
		}
		mat_out[k] = temp;
	}
}

// T should be generic C++ native types
template <typename T, int nrows, int ncols>
void arma2dlib(const arma::Cube<T> &mat_in, std::vector<dlib::matrix<T, nrows, ncols>> &mat_out)
{
	//int nr = mat_in.n_rows;
	//int nc = mat_in.n_cols;
	int nr = nrows;
	int nc = ncols;
	int nchannels = mat_in.n_slices;

	mat_out.resize(nchannels);

	dlib::matrix<T, nrows, ncols> temp;

	for (int k = 0; k < nchannels; k++)
	{
		for (int j = 0; j < nc; j++)
		{
			for (int i = 0; i < nr; i++)
				temp(i, j) = mat_in.at(i,j,k);
		}
		mat_out[k] = temp;
	}
}

// RGB color image
// This function makes a deep copy.
void arma2dlib(const arma::Cube<unsigned char> &mat_in, dlib::array2d<dlib::rgb_pixel> &mat_out)
{

	int nrows = mat_in.n_rows;
	int ncols = mat_in.n_cols;
	const int nchannels = 3;

	mat_out.set_size(nrows, ncols);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].red = mat_in.at(i,j,0);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].green = mat_in.at(i,j,1);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].blue = mat_in.at(i,j,2);
}

// BGR color image
// This function makes a deep copy.
void arma2dlib(const arma::Cube<unsigned char> &mat_in, dlib::array2d<dlib::bgr_pixel> &mat_out)
{

	int nrows = mat_in.n_rows;
	int ncols = mat_in.n_cols;
	const int nchannels = 3;

	mat_out.set_size(nrows, ncols);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].blue= mat_in.at(i,j,0);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].green = mat_in.at(i,j,1);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].red = mat_in.at(i,j,2);
}

// HSI color image
// This function makes a deep copy.
void arma2dlib(const arma::Cube<unsigned char> &mat_in, dlib::array2d<dlib::hsi_pixel> &mat_out)
{

	int nrows = mat_in.n_rows;
	int ncols = mat_in.n_cols;
	const int nchannels = 3;

	mat_out.set_size(nrows, ncols);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].h = mat_in.at(i,j,0);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].s = mat_in.at(i,j,1);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].i = mat_in.at(i,j,2);
}

// LAB color image
// This function makes a deep copy.
void arma2dlib(const arma::Cube<unsigned char> &mat_in, dlib::array2d<dlib::lab_pixel> &mat_out)
{

	int nrows = mat_in.n_rows;
	int ncols = mat_in.n_cols;
	const int nchannels = 3;

	mat_out.set_size(nrows, ncols);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].l = mat_in.at(i,j,0);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].a = mat_in.at(i,j,1);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].b = mat_in.at(i,j,2);
}

// RGBA color image
// This function makes a deep copy.
void arma2dlib(const arma::Cube<unsigned char> &mat_in, dlib::array2d<dlib::rgb_alpha_pixel> &mat_out)
{
	int nrows = mat_in.n_rows;
	int ncols = mat_in.n_cols;
	const int nchannels = 4;

	mat_out.set_size(nrows, ncols);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].red = mat_in.at(i,j,0);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].green = mat_in.at(i,j,1);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].blue = mat_in.at(i,j,2);

	for (int j = 0; j < ncols; j++)
		for (int i = 0; i < nrows; i++)
			mat_out[i][j].alpha = mat_in.at(i,j,3);
}





#endif