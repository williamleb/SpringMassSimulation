/**
 * @file ParticleSystem.cpp
 *
 * @brief Système de particules de type masse-ressort
 *
 * Nom: William Lebel
 * Email : william.lebel.1@ens.etsmtl.ca
 *
 */

#include "ParticleSystem.h"

using namespace gti320;

const double ParticleSystem::gravitationalConstant = 9.81;

/**
 * Calcule des forces qui affectent chacune des particules.
 *
 * Étant donné une particule p, la force est stockée dans le vecteur p.f
 * Les forces prisent en compte sont : la gravité et la force des ressorts.
 */
void ParticleSystem::computeForces()
{
	// Calcul de la force gravitationnelle sur chacune des particules
	for (Particle& particle : m_particles)
	{
		particle.f.x() = 0.0;
		particle.f.y() = -gravitationalConstant * particle.m;
	}

	// Calcul de la force de chaque ressort pour ses particules associées
	for (const Spring& spring : m_springs)
	{
		auto& firstParticle = m_particles[spring.index0];
		auto& secondParticle = m_particles[spring.index1];

		auto differenceVector = secondParticle.x - firstParticle.x;
		auto distance = differenceVector.norm();
		auto direction = differenceVector * (1.0 / distance);

		auto force = spring.k * (distance - spring.l0);
		auto forceVector = direction * force;

		firstParticle.f = firstParticle.f + forceVector;
		secondParticle.f = secondParticle.f - forceVector;
	}
}

/**
 * Assemble les données du système dans les vecteurs trois vecteurs d'état outPos,
 * outVel et outForce.
 */
void ParticleSystem::pack(Vector<double, Dynamic>& outPos,
                          Vector<double, Dynamic>& outVel,
                          Vector<double, Dynamic>& outForce)
{
	const int numberOfParticles = static_cast<int>(m_particles.size());
	outPos.resize(numberOfParticles * 2);
	outVel.resize(numberOfParticles * 2);
	outForce.resize(numberOfParticles * 2);

	for (auto i = 0; i < numberOfParticles; ++i)
	{
		auto particle = m_particles[i];

		outPos(2 * i) = particle.x.x();
		outPos(2 * i + 1) = particle.x.y();

		outVel(2 * i) = particle.v.x();
		outVel(2 * i + 1) = particle.v.y();

		outForce(2 * i) = particle.f.x();
		outForce(2 * i + 1) = particle.f.y();
	}
}

/**
 * Copie les données des vecteurs d'états dans le particules du système.
 */
void ParticleSystem::unpack(const Vector<double, Dynamic>& pos,
                            const Vector<double, Dynamic>& vel)
{
	ASSERTF(pos.size() == static_cast<int>(m_particles.size()) * 2, "Trying to unpack positions with a particle system of uncompatible size (%d vs %d)", pos.size(), static_cast<int>(m_particles.size()));
	ASSERTF(vel.size() == static_cast<int>(m_particles.size()) * 2, "Trying to unpack velocities with a particle system of uncompatible size (%d vs %d)", vel.size(), static_cast<int>(m_particles.size()));

	for (auto i = 0; i < static_cast<int>(m_particles.size()); ++i)
	{
		auto& particle = m_particles[i];

		particle.x.x() = pos(2 * i);
		particle.x.y() = pos(2 * i + 1);

		particle.v.x() = vel(2 * i);
		particle.v.y() = vel(2 * i + 1);
	}
}


/**
 * Construction de la matrice de masses.
 */
void ParticleSystem::buildMassMatrix(Matrix<double, Dynamic, Dynamic>& outMassMatrix)
{
	const int numberOfParticles = static_cast<int>(m_particles.size());
	const int dimensions = 2 * numberOfParticles;
	outMassMatrix.resize(dimensions, dimensions);
	outMassMatrix.setZero();

	for (int i = 0; i < numberOfParticles; ++i)
	{
		auto particle = m_particles[i];
		auto particleMass = !particle.fixed ? particle.m : std::numeric_limits<double>::max();

		outMassMatrix(2 * i, 2 * i) = particleMass;
		outMassMatrix(2 * i + 1, 2 * i + 1) = particleMass;
	}
}


/**
 * Construction de la matrice de rigidité.
 */
void ParticleSystem::buildDfDx(Matrix<double, Dynamic, Dynamic>& outDfDxMatrix)
{
	const int numberOfParticles = m_particles.size();
	const int dimensions = 2 * numberOfParticles;
	outDfDxMatrix.resize(dimensions, dimensions);
	outDfDxMatrix.setZero();

	// Pour chaque ressort...
	for (const Spring& spring : m_springs)
	{
		auto firstParticle = m_particles[spring.index0];
		auto secondParticle = m_particles[spring.index1];

		auto differenceVector = secondParticle.x - firstParticle.x;
		auto squaredDistance = differenceVector.squaredNorm();
		auto distance = sqrt(squaredDistance);

		auto springContribution = this->dyadicProduct(differenceVector, differenceVector);

		auto alphaCoefficient = spring.k * (1.0 - (spring.l0 / distance));
		auto dyadicProductCoefficient = spring.k * (spring.l0 / (distance * squaredDistance));

		springContribution(0, 0) = springContribution(0, 0) * dyadicProductCoefficient + alphaCoefficient;
		springContribution(1, 1) = springContribution(1, 1) * dyadicProductCoefficient + alphaCoefficient;
		springContribution(0, 1) = springContribution(0, 1) * dyadicProductCoefficient;
		springContribution(1, 0) = springContribution(1, 0) * dyadicProductCoefficient;

		// On ajoute la contribution à la diagonale
		outDfDxMatrix.block(spring.index0 * 2, spring.index0 * 2, 2, 2) -= springContribution;
		outDfDxMatrix.block(spring.index1 * 2, spring.index1 * 2, 2, 2) -= springContribution;

		// On ajoute la contribution négative aux autres endroits de la matrice
		outDfDxMatrix.block(spring.index0 * 2, spring.index1 * 2, 2, 2) = springContribution;
		outDfDxMatrix.block(spring.index1 * 2, spring.index0 * 2, 2, 2) = springContribution;
	}
}

Matrix<double, 2, 2> ParticleSystem::dyadicProduct(const Vector2d& left, const Vector2d& right) const
{
	Matrix<double, 2, 2> dyadicProduct;

	dyadicProduct(0, 0) = left.x() * right.x();
	dyadicProduct(0, 1) = left.x() * right.y();

	dyadicProduct(1, 0) = left.y() * right.x();
	dyadicProduct(1, 1) = left.y() * right.y();

	return dyadicProduct;
}
