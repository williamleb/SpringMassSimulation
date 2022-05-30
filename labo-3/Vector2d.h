#pragma once

/**
 * @file Vector2d.h
 *
 * @brief Vecteur 2D
 *
 * Nom: William Lebel
 * Email : william.lebel.1@ens.etsmtl.ca
 *
 */

#include "Math3D.h"

namespace gti320
{
	/**
	 * Spécialisation de template pour le les vecteurs 2D
	 */
	template <typename Scalar>
	class Vector<Scalar, 2> : public MatrixBase<Scalar, 2, 1>
	{
	public:
		/**
		 * Constructeur par défaut
		 */
		Vector() : MatrixBase<Scalar, 2, 1>()
		{
		}

		/**
		 * Constructeur à partir de la dimension. Afin d'assurer la compatibilité
		 * avec la version non-spécialisée de la classe.
		 *
		 * @param size doit être 2
		 */
		Vector(int size) : MatrixBase<Scalar, 2, 1>()
		{
			ASSERT(size == 2, "Vector 2's size must be greater than 0");
		}

		/**
		 * Constructeur à partir des coordonnées x et y
		 */
		Vector(Scalar x, Scalar y) : MatrixBase<Scalar, 2, 1>()
		{
			(*this)(0) = x;
			(*this)(1) = y;
		}

		/**
		 * Constructeur de copie
		 */
		Vector(const Vector& other) : MatrixBase<Scalar, 2, 1>(other)
		{
		}

		/**
		 * Destructeur
		 */
		~Vector()
		{
		}

		/**
		 * Opérateur de copie
		 */
		Vector& operator=(const Vector& other)
		{
			this->m_storage = other.m_storage;
			return *this;
		}

		/**
		 * Accesseur à une entrée du vecteur (lecture seule)
		 */
		Scalar operator()(int i) const
		{
			ASSERT(0 <= i && i < 2, "Trying to access a value in a vector2 that is not 0 or 1");
			return *(this->m_storage.data() + i);
		}

		/**
		 * Accesseur à une entrée du vecteur (lecture et écriture)
		 */
		Scalar& operator()(int i)
		{
			ASSERT(0 <= i && i < 2, "Trying to access a value in a vector2 that is not 0 or 1");
			return *(this->m_storage.data() + i);
		}


		/**
		 * Produit scalaire de *this et other
		 */
		inline Scalar dot(const Vector& other) const
		{
			return (*this)(0) * other(0) + (*this)(1) * other(1);
		}

		/**
		 * Retourne la norme euclidienne à la 2 du vecteur
		 */
		inline Scalar squaredNorm() const
		{
			return dot(*this);
		}

		/**
		 * Norme eulidienne du vecteur
		 */
		inline Scalar norm() const
		{
			return sqrt(dot(*this));
		}

		/**
		 * Retourne le vecteur normalisé
		 */
		inline Scalar normalize() const
		{
			return (*this) * (static_cast<Scalar>(1) / norm());
		}

		/**
		 * Retourne la distance euclidienne du vecteur avec un autre vecteur
		 */
		inline Scalar squaredDistance(const Vector& other) const
		{
			ASSERT(this->size() == other.size(), "Trying to find the distance between two vectors that are not the same size");

			Scalar distanceX = other.x() - this->x();
			Scalar distanceY = other.y() - this->y();

			return distanceX * distanceX + distanceY * distanceY;
		}

		/**
		 * Retourne la distance euclidienne du vecteur avec un autre vecteur
		 */
		inline Scalar distance(const Vector& other) const
		{
			return sqrt(this->squaredDistance(other));
		}

		inline Scalar x() const { return (*this)(0); }
		inline Scalar& x() { return (*this)(0); }
		inline Scalar y() const { return (*this)(1); }
		inline Scalar& y() { return (*this)(1); }
	};

	typedef Vector<double, 2> Vector2d;
}
