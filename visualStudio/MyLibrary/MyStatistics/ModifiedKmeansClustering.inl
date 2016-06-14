namespace MyStatistics{

	inline ModifiedKmeansClustering::ModifiedKmeansClustering(const size_t cluster_num_, const std::vector<double>& radiusMinimum)
		: cluster_num(cluster_num_)
	{
		clusters = std::vector<ModifiedKmeansClusters>(0);
		for (size_t c = 0; c < cluster_num; ++c)
			clusters.push_back(ModifiedKmeansClusters(c, radiusMinimum));
		
		data = std::vector<KmeansData>(0);
		rnd = gsl_rng_alloc(gsl_rng_default);
	}
	inline ModifiedKmeansClustering::~ModifiedKmeansClustering(void){
		gsl_rng_free(rnd);
	}

	inline double ModifiedKmeansClustering::runMultiple(){
		AICmin = run();
		for (size_t i = 0; i < data.size(); ++i)
			data[i].setClusterIndexRslt();

		for (size_t t = 0; t < 5; ++t){
			double AIC = run();
			if (AIC < AICmin){
				AICmin = AIC;
				for (size_t i = 0; i < data.size(); ++i)
					data[i].setClusterIndexRslt();
			}
		}
/*
#ifdef DEBUG
		std::vector<double> avg(data.size());
		std::vector<double> sigma(data.size());
		std::vector<double> cluster(data.size());
		for (size_t i = 0; i < data.size(); ++i){
			avg[i] = data[i].avg;
			sigma[i] = data[i].sigma;
			cluster[i] = 1.0*data[i].clusterIndexRslt;
		}
		IGORdata::write_itx(avg, sigma, "avg_sigma.itx", "avg", "sigma");
		IGORdata::write_itx(cluster, "cluster.itx", "cluster");
#endif
*/
		return AICmin;
	}


	inline double ModifiedKmeansClustering::run(){
		initialize();

		for (size_t i = 0; i < 100 &&
			updateClusters();
			++i){
		}

		double AIC = 0.0;
		for (size_t i = 0; i < data.size(); ++i){
			size_t c = data[i].clusterIndex;
			AIC += clusters[c].getLiklihood(data[i]);
		}
		AIC = -2.0*log(AIC) - 2.0*2.0;
		return AIC;
	}


	inline bool ModifiedKmeansClustering::initialize(){

		if (data.size() == 0) return false;

		for (size_t i = 0; i < data.size(); ++i)
			data[i].clusterIndex = 0;

		clusters[0].set(data);

		std::vector<double> initialAvg(clusters[0].numDims);
		for (size_t d = 0; d < clusters[0].numDims; ++d)
			initialAvg[d] = gsl_ran_gaussian(rnd, sqrt(clusters[0].getSigma(d))) 
							          + gsl_vector_get(clusters[0].avg, d);
		
		clusters[0].setInitialAvg(initialAvg);


		//	initial guess of the other clusters
		std::vector<double> squaredRootWeight(data.size());
		double total = 0.0;
		for (size_t i = 0; i < data.size(); ++i){
			squaredRootWeight[i] = total;
			total += clusters[0].getDistanceSq(data[i]);
		}
		for (size_t c = 1; c < cluster_num; ++c){
			double value = gsl_rng_uniform(rnd) * total;
			for (size_t i = 0; i < data.size(); ++i){
				if (squaredRootWeight[i] > value){
					clusters[c].setInitialAvg(data[i].data);
					//	copy cov matrix
					gsl_matrix_memcpy(clusters[c].cov, clusters[0].cov);
					clusters[c].detCovSqrtInv = clusters[0].detCovSqrtInv;
					break;
				}
			}
		}
		return true;
	}
	
	
	inline bool ModifiedKmeansClustering::updateBelongingClusters(){
		bool isUpdated = false;
		for (size_t i = 0; i < data.size(); ++i){
			size_t cindex = data[i].clusterIndex;
			double liklihood = clusters[cindex].getLiklihood(data[i]);

			for (size_t c = 0; c < cluster_num; ++c){
				double liklihoodNew = clusters[c].getLiklihood(data[i]);
				if (liklihoodNew > liklihood && c != cindex){
					liklihood = liklihoodNew;
					data[i].clusterIndex = c;
					isUpdated = true;
				}
			}
		}
		return isUpdated;
	}


	inline bool ModifiedKmeansClustering::updateClusters(){
		for (size_t c = 0; c < cluster_num; ++c){
			if (!updateBelongingClusters() || !clusters[c].set(data))
				return false;
		}
		return true;
	}

	//-----------------------------------------
	//
	//			ModifiedKmeansClusters
	//
	//-----------------------------------------
	inline ModifiedKmeansClusters::ModifiedKmeansClusters(void)	//	unused constructor
		:KmeansClusters()
	{}//	unused constructor

	inline ModifiedKmeansClusters::ModifiedKmeansClusters
		(const size_t clusterIndex_, const std::vector<double>& radiusMinimum_)
		:
		KmeansClusters(clusterIndex_, radiusMinimum_)
	{	}

	inline ModifiedKmeansClusters::~ModifiedKmeansClusters(void)
	{	}

	inline ModifiedKmeansClusters::ModifiedKmeansClusters(const ModifiedKmeansClusters& obj)
		:
		KmeansClusters(obj)
	{	}

	inline ModifiedKmeansClusters& ModifiedKmeansClusters::operator = (const ModifiedKmeansClusters& obj)
	{
		copy(obj);
		return *this;
	}

	//	return false if covariant matrix cannot be calculated
	inline bool ModifiedKmeansClusters::set(const std::vector<KmeansData>& data)
	{
		gsl_vector_set_zero(avg);
		size_t count = 0;
		//	evaluate the average value of the data
		for (size_t i = 0; i < data.size(); ++i){
			if (data[i].clusterIndex == clusterIndex){
				for (size_t d = 0; d < numDims; ++d){
					gsl_vector_set(avg, d, gsl_vector_get(avg, d) + data[i].data[d]);
					count++;
				}
			}
		}

		if (count == 0) return false;
		for (size_t d = 0; d < numDims; ++d)
			gsl_vector_set(avg, d, gsl_vector_get(avg, d) / count);

		//	evaluate the covariant matrix
		for (size_t i = 0; i < data.size(); ++i){
			if (data[i].clusterIndex == clusterIndex){
				for (size_t d0 = 0; d0 < numDims; ++d0){
					double tmp0 = data[i].data[d0] - gsl_vector_get(avg, d0);
					for (size_t d1 = 0; d1 < numDims; ++d1){
						double tmp1 = data[i].data[d1] - gsl_vector_get(avg, d1);
						gsl_matrix_set(cov, d0, d1, tmp0*tmp1 + gsl_matrix_get(cov, d0, d1));
					}
				}
			}
		}

		for (size_t d0 = 0; d0 < numDims; ++d0){
			for (size_t d1 = 0; d1 < numDims; ++d1){
				if (d0 == d1)		//	leveling of the minimum 
					gsl_matrix_set(cov, d0, d0, gsl_matrix_get(cov, d0, d1) / count + radiusMinimum[d0] * radiusMinimum[d0]);
				else
					gsl_matrix_set(cov, d0, d1, gsl_matrix_get(cov, d0, d1) / count);
			}
		}

		//	invert the matrix cov using the LU decomposition
		gsl_matrix_memcpy(LU, cov);
		int signum = 0;

		gsl_linalg_LU_decomp(LU, perm, &signum);
		double det = gsl_linalg_LU_det(LU, signum);
		detCovSqrtInv = 1.0 / sqrt(det);
		gsl_linalg_LU_invert(LU, perm, cov);
		
		return true;
	}

};
