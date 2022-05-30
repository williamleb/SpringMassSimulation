#pragma once

/**
 * @file Solvers.hpp
 *
 * @brief Implémentation de plusieurs algorihtmes de solveurs pour un système
 *        d'équations linéaires
 *
 * Nom: William Lebel
 * Email : william.lebel.1@ens.etsmtl.ca
 *
 */

#include <omp.h>

#include "Math3D.h"
#include <stdio.h>

namespace gti320
{
	// Identification des solveurs
	enum eSolverType { kNone, kJacobi, kGaussSeidel, kCholesky };

	// Paramètres de convergences pour les algorithmes itératifs
	static const double epsilon = 1e-4;
	static const double tau = 1e-5;

	/**
	 * Résout Ax = b avec la méthode de Jacobi
	 */
	static void jacobi(const Matrix<double, Dynamic, Dynamic>& A,
	                   const Vector<double, Dynamic>& b,
	                   Vector<double, Dynamic>& outSolution, int k_max)
	{
		ASSERT(A.rows() == A.cols(), "Trying to apply Jacobi solver with a non square matrix");
		ASSERT(b.size() == A.rows(), "Trying to apply Jacobi solver with a vector of size incompatible with the matrix");

		outSolution = b;

		// z: partialSolution
		auto partialSolution = b;

		auto size = A.rows();

		auto numberOfIterations = 0;
		Vector<double, Dynamic> lastSolution;
		do
		{
			lastSolution = outSolution;

			#pragma omp parallel for
			for (auto i = 0; i < size; ++i)
			{
				auto partialSolutionElement = b(i);

				for (auto j = 0; j < i; ++j)
				{
					partialSolutionElement -= A(i, j) * outSolution(j);
				}
				for (auto j = i + 1; j < size; ++j)
				{
					partialSolutionElement -= A(i, j) * outSolution(j);
				}

				partialSolutionElement /= A(i, i);
				partialSolution(i) = partialSolutionElement;
			}

			outSolution = partialSolution;
			++numberOfIterations;

		} while (numberOfIterations < k_max && (outSolution - lastSolution).norm() / outSolution.norm() > tau && (A * outSolution - b).norm() / b.norm() > epsilon);
	}


	/**
	 * Résout Ax = b avec la méthode Gauss-Seidel
	 */
	static void gaussSeidel(const Matrix<double, Dynamic, Dynamic>& A,
	                        const Vector<double, Dynamic>& b,
	                        Vector<double, Dynamic>& outSolution, int k_max)
	{
		ASSERT(A.rows() == A.cols(), "Trying to apply Gauss-Seidel solver with a non square matrix");
		ASSERT(b.size() == A.rows(), "Trying to apply Gauss-Seidel solver with a vector of size incompatible with the matrix");

		outSolution = b;

		auto size = A.rows();
		auto numberOfIterations = 0;
		Vector<double, Dynamic> lastSolution;
		do
		{
			lastSolution = outSolution;

			for (auto i = 0; i < size; ++i)
			{
				outSolution(i) = b(i);

				for (auto j = 0; j < i; ++j)
				{
					outSolution(i) -= A(i, j) * outSolution(j);
				}
				for (auto j = i + 1; j < size; ++j)
				{
					outSolution(i) -= A(i, j) * outSolution(j);
				}

				outSolution(i) /= A(i, i);
			}

			++numberOfIterations;

		} while (numberOfIterations < k_max && (outSolution - lastSolution).norm() / outSolution.norm() > tau && (A * outSolution - b).norm() / b.norm() > epsilon);
	}

	/**
	 * Résout Ax = b avec la méthode de Cholesky
	 * @param A A
	 * @param b b
	 * @param outSolution x
	 */
	static void cholesky(const Matrix<double, Dynamic, Dynamic>& A,
	                     const Vector<double, Dynamic>& b,
	                     Vector<double, Dynamic>& outSolution)
	{
		ASSERT(A.rows() == A.cols(), "Trying to apply Cholesky factorization to a non square matrix");
		ASSERT(b.size() == A.rows(), "Trying to apply Cholesky solver with a vector of size incompatible with the matrix");

		outSolution.resize(b.size());

		auto size = A.rows();

		// Calcul de la matrice L de la factorisation de Cholesky
		// L: cholesky
		Matrix<double, Dynamic, Dynamic> cholesky(A.rows(), A.cols());
		cholesky.setZero();

		for (auto i = 0; i < size; ++i)
		{
			for (auto k = 0; k <= i; ++k)
			{
				auto sum = 0.0;
				for (auto j = 0; j < k; ++j)
				{
					sum += cholesky(i, j) * cholesky(k, j);
				}

				if (i == k)
				{
					cholesky(i, k) = sqrt(A(i, i) - sum);
				}
				else
				{
					cholesky(i, k) = (A(i, k) - sum) / cholesky(k, k);
				}
			}
		}

		// Résout Ly = b
		// y: partialSolution
		Vector<double, Dynamic> partialSolution(b.size());
		for (auto i = 0; i < size; ++i)
		{
			partialSolution(i) = b(i);

			for (auto j = 0; j < i; ++j)
			{
				partialSolution(i) -= cholesky(i, j) * partialSolution(j);
			}

			partialSolution(i) /= cholesky(i, i);
		}

		// Résout L^t x = y
		for (auto i = size - 1; i >= 0; --i)
		{
			outSolution(i) = partialSolution(i);

			for (auto j = i + 1; j < size; ++j)
			{
				outSolution(i) -= cholesky(j, i) * outSolution(j);
			}

			outSolution(i) /= cholesky(i, i);
		}
	}
}
