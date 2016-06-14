namespace MeanShiftClustering{
	//-----------------------------------------------
	//		static functions of Cluseterize
	//-----------------------------------------------
	template < class Cor_ >
	inline double Clusterize< Cor_ >::kernel(double deltap_squared){
		return exp(-0.5*deltap_squared);
	}
	template < class Cor_ >
	inline double Clusterize< Cor_ >::dkernel_dx(double deltap_squared){
		return exp(-0.5*deltap_squared)*0.5;
	}

	//-----------------------------------------------
	//		member functions of Clusterize
	//-----------------------------------------------

	template < class Cor_ >
	inline void Clusterize< Cor_ >::setData(const std::vector<ClusteringData <Cor_> >& data_){
		allocate(data_.size());
		for (size_t i = 0; i < size; ++i)
			data[i] = data_[i];
	}
	template < class Cor_ >
	inline void Clusterize< Cor_ >::setData(const ClusteringData< Cor_ >* data_, const size_t size_){
		allocate(size_);
		for (size_t i = 0; i < size; ++i)
			data[i] = data_[i];
	}
	template < class Cor_ >
	inline void Clusterize< Cor_ >::setOneData(const size_t dataIndex, const ClusteringData< Cor_ >& data_){
		data[dataIndex] = data_;
	}


	template < class Cor_ >
	inline void Clusterize< Cor_ >::setThreshold(const double threshold_iteration_, const double threshold_clusterize_){
		threshold_iteration = threshold_iteration_;
		threshold_clusterize = threshold_clusterize_;
	}

	//	core function
	template < class Cor_ >
	inline void Clusterize<Cor_>::clusterize(){
		for (size_t i = 0; i < size; ++i){
			if (data[i].isValid)
				meanShiftOne(data[i]);
		}
		labelling();
	}

	template < class Cor_ >
	inline void Clusterize< Cor_ >::meanShiftOne(ClusteringData< Cor_ >& obj){
		double threshold_squared = threshold_iteration*threshold_iteration;

		Cor_ p_now(obj.x);
		Cor_ p_prev(obj.x);
		do{
			p_prev = p_now;
			p_now = updateMean(p_prev);
		} while (p_now.getDeltaPSquared(p_prev) > threshold_squared);
		obj.xc = p_now;
		return;
	}

	template < class Cor_ >
	inline Cor_ Clusterize< Cor_ >::updateMean(Cor_& p_start)const{
		Cor_ sigma_gi_xi;
		sigma_gi_xi.setZeros();
		double sigma_gi = 0.0;

		for (size_t i = 0; i < size; ++i){
			if (data[i].isValid){
				double gi = dkernel_dx(p_start.getDeltaPSquared(data[i].x));
				sigma_gi += gi;
				sigma_gi_xi.add(data[i].x*gi);
			}
		}
		return sigma_gi_xi / sigma_gi;
	}

	
	template < class Cor_ >
	inline void Clusterize< Cor_ >::labelling(){
		std::vector<std::vector<size_t>> clusterSize = std::vector<std::vector<size_t>>(0);

		int label_now = -1;
		for (size_t i = 0; i < size; ++i){
			if (data[i].isValid && data[i].label < 0){
				++label_now;
				data[i].label = label_now;
				std::vector<size_t> cluster_index_one(0);
				cluster_index_one.push_back(i);
				for (size_t j = i + 1; j < size; ++j){
					if (data[j].isValid){
						if (data[j].label < 0 && data[i].xc.getDeltaPSquared(data[j].xc) < threshold_clusterize){
							data[j].label = label_now;
							cluster_index_one.push_back(j);
						}
					}
				}
				clusterSize.push_back(cluster_index_one);
			}
		}
		clusters.allocate(label_now+1);
		for (size_t i = 0; i < label_now + 1; ++i){
			clusters.memberClusters[i].allocate(clusterSize[i].size());
			for (size_t j = 0; j < clusters.memberClusters[i].Size(); ++j){
				clusters.memberClusters[i].memberData[j] = &(data[i]);
			}
		}
	}

	//---	memory management	---
	template < class Cor_ >
	inline void Clusterize< Cor_ >::allocate(const size_t size_){
		if (isAllocated)	free();
		size = size_;
		data = new ClusteringData< Cor_ >[size];
		isAllocated = true;
	}
	template < class Cor_ >
	inline void Clusterize< Cor_ >::free(){
		if (isAllocated) delete[] data;
		clusters.free();
	}


	//-----------------------------------------------
	//		member functions of Cluster
	//-----------------------------------------------
	//---	memory management	---
	template < class Cor_ >
	inline void Cluster< Cor_ >::allocate(const size_t size_){
		if (isAllocated)	free();
		size = size_;
		memberData = new ClusteringData< Cor_ >*[size];
		isAllocated = true;
	}
	template < class Cor_ >
	inline void Cluster< Cor_ >::free(){
		if (isAllocated) delete[] memberData;
	}

	//-----------------------------------------------
	//		member functions of Clusters
	//-----------------------------------------------
	//---	memory management	---
	template < class Cor_ >
	inline void Clusters< Cor_ >::allocate(const size_t size_){
		if (isAllocated)	free();
		size = size_;
		memberClusters = new Cluster< Cor_ >[size];
		isAllocated = true;
	}
	template < class Cor_ >
	inline void Clusters< Cor_ >::free(){
		for (size_t i = 0; i < size; ++i)
			memberClusters[i].free();
		if (isAllocated) delete[] memberClusters;
	}

};