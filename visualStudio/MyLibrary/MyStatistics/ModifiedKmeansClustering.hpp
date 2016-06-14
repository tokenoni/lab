#ifndef __MODIFIED_ModifiedKmeans_CLUSTERING_HPP__
#define __MODIFIED_ModifiedKmeans_CLUSTERING_HPP__

#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "KmeansClustering.hpp"
#define GSL_HAVE_INLINE

namespace MyStatistics{
	class ModifiedKmeansClusters;


	class ModifiedKmeansClustering : public KmeansClustering {
	public:
		ModifiedKmeansClustering(const size_t cluster_num_, const std::vector<double>& radiusMinimum_);
		~ModifiedKmeansClustering();

		//	returns the AIC value
		double runMultiple();

	private:
		//	returns the AIC value
		double run();
		bool initialize();
		const size_t cluster_num;

		bool updateBelongingClusters();
		bool updateClusters();

		std::vector<KmeansData> data;
		std::vector<ModifiedKmeansClusters> clusters;
		gsl_rng* rnd;

		double AICmin;

	};

	class ModifiedKmeansClusters : public KmeansClusters{
	public:
		//	clusterIndex:	id of each object
		//	radiusMinimum_:	minimum of the radius for each dimension. 
		ModifiedKmeansClusters(void);
		ModifiedKmeansClusters
			(const size_t clusterIndex, 
			 const std::vector<double>& radiusMinimum_);
		~ModifiedKmeansClusters();
		ModifiedKmeansClusters(const ModifiedKmeansClusters& obj);
		ModifiedKmeansClusters& operator = (const ModifiedKmeansClusters& obj);
		
		//	return false if covariant matrix cannot be calculated
		bool set(const std::vector<KmeansData>& data);

	};


};

#include "ModifiedKmeansClustering.inl"
#endif