#ifndef __MEANSHIFT_CLUSTERING_HPP__
#define __MEANSHIFT_CLUSTERING_HPP__

#include <vector>

namespace MeanShiftClustering{
	//	template class must have the following members;
	class Coordinate{
	public:
		virtual void setZeros();
		virtual double getDeltaPSquared(const Coordinate src)const;
		virtual Coordinate add(const Coordinate src);
		virtual Coordinate operator* (const Coordinate src)const;
		virtual Coordinate operator* (const double src)const;
		virtual Coordinate operator/ (const double src)const;
		virtual Coordinate& operator= (const Coordinate src);
	};

	template < class Cor_ > class Clusterize;
	template < class Cor_ > class ClusteringData;
	template < class Cor_ > class Clusters;
	template < class Cor_ > class Cluster;


	template < class Cor_ >
	class ClusteringData{
	public:
		ClusteringData(void){};
		ClusteringData(Cor_& x_){
			x = x_; isValid = true; label = -1;
		};
	
		Cor_ x;
		int label;
		bool isValid;	//	false if y<0 

		Cor_ xc;	//	iterated point;
	};

	template < class Cor_ >
	class Clusterize{
	public:
		Clusterize(void){ isAllocated = false; }

		//	1. set data
		//	option 1:	data is given as vector
		void setData(const std::vector<ClusteringData<Cor_>>& data_);
		//	option 2:	data is given as pointer list
		void setData(const ClusteringData<Cor_>* data_, const size_t size_);
		//	option 3:	data is given separately
		void setSize(const size_t size_){ allocate(size_); };
		void setOneData(const size_t dataIndex, const ClusteringData<Cor_>& data_);


		//	2. set the threshold.
		void setThreshold(const double threshold_iteration = 0.01, const double threshold_clusterize = 1.0);
	 
		//	3. clusterize 
		void clusterize();

		//	tip functions
		size_t Size()const{ return size; }

		//	access to data
		const ClusteringData< Cor_ >& getData(const size_t dataIndex)const{
			return data[dataIndex];
		}

		//	access to results
		const ClusteringData< Cor_ >& getClusterizedData(const size_t clusterIndex, const size_t memberIndex)const{
			return clusters.memberClusters[clusterIndex].memberData[memberIndex];
		}
		const Cluster< Cor_ >& getCluster(const size_t clusterIndex)const{
			return clusters.memberClusters[clusterIndex];
		}
		const Cluster< Cor_ >& getClusterFromData(const size_t dataIndex)const{
			if (data[dataIndex].isValid)
				return clusters.memberClusters[data[dataIndex].label];
			else {
				Cluster< Cor_ > tmp;
				return tmp;
			}
		}

	protected:
		ClusteringData< Cor_ >* data;
		Clusters< Cor_ > clusters;

		size_t size;
		double threshold_iteration, threshold_clusterize;

		void meanShiftOne(ClusteringData< Cor_>& obj);
		Cor_ updateMean(Cor_& p_start)const;
		Cor_ getGradF(double x, double y)const;
		void labelling();
		void save();

		//	for memory management
		void allocate(const size_t size_);
		void free();
		bool isAllocated;
		
	public:
		//	these functions can be overriden
		static double kernel(double deltap_squared);
		static double dkernel_dx(double deltap_squared);

	};

	template < class Cor_ >
	class Clusters{
		friend ClusteringData<Cor_>;
		friend Clusterize<Cor_>;
	public:
		Clusters(void){ isAllocated = false; }
		size_t Size()const{ return size; }
	private:
		void allocate(const size_t size_);
		
		Cluster< Cor_ >* memberClusters;
		size_t size;
		bool isAllocated;
		void free();
	};

	template < class Cor_ >
	class Cluster{
		friend Clusterize<Cor_>;
		friend Clusters<Cor_>;
	public:
		Cluster(void){ isAllocated = false; }
		size_t Size()const{ return size; }
	protected:
		void allocate(const size_t size_);
		size_t size;
		ClusteringData<Cor_>** memberData;	//	pointer to the original data. Original data is always stored in the "Clusterize" class
		bool isAllocated;
		void free();
	};
};

#include "MeanShiftClustering.inl"
#endif