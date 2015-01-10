#include<iostream>
#include<vector>

#include<cmath>
#include<cstdlib>
#include<cfloat>

#include<Eigen/Dense>

#include<pcl/common/common_headers.h>
#include<pcl/common/centroid.h>
#include<pcl/common/distances.h>
#include<pcl/visualization/pcl_visualizer.h>

#include<boost/thread/thread.hpp>
#include<boost/tuple/tuple.hpp>

template<class PointT>
class MLSDeform{
private:
	const typename pcl::PointCloud<PointT>::Ptr cloud;
	const typename pcl::PointCloud<PointT>::Ptr original;
	const typename pcl::PointCloud<PointT>::Ptr shifted;
	const float alpha;

	const float eucdist(const pcl::PointXYZ &a,const pcl::PointXYZ &b){
		return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
	}

	inline PointT E2P(const Eigen::MatrixXf &in){
		return PointT(in(0,0),in(1,0),in(2,0));
	}

	inline Eigen::Matrix<float,3,1> P2E(const pcl::PointXYZ &in){
		Eigen::Matrix<float,3,1> out;
		out<<in.x,in.y,in.z;
		return out;
	}

	Eigen::MatrixXf PointCloud2Eigen(const typename pcl::PointCloud<PointT>::Ptr cloud){
		Eigen::MatrixXf out(3,cloud->size());
		int j=0;
		for(typename pcl::PointCloud<PointT>::iterator i=cloud->begin();i<cloud->end();i++){
			Eigen::MatrixXf temp(3,1);
			temp<<i->x,i->y,i->z;
			out.col(j++) = temp;
		}
		return out;
	}

public:
	MLSDeform(const typename pcl::PointCloud<PointT>::Ptr cloud_,
	                               const typename pcl::PointCloud<PointT>::Ptr original_,
																 const typename pcl::PointCloud<PointT>::Ptr shifted_,
																 const float alpha)
																:cloud(cloud_),
																 original(original_),
																 shifted(shifted_),
																 alpha(alpha){}

	typename pcl::PointCloud<PointT>::Ptr deform(){
		const int cloudSize = cloud->size();
		const long controlSize = original->size();
		typename pcl::PointCloud<PointT>::Ptr final(new pcl::PointCloud<PointT>);
		final->reserve(cloudSize);
		Eigen::MatrixXf weights(cloudSize,controlSize);

		for(int i=0;i<cloudSize;i++){//Initialise the weight matrix
			for(int j=0;j<controlSize;j++){
				const float distance = eucdist(cloud->points[i],original->points[j]);
				if(distance>0.001){
					weights(i,j) = pow(distance,-2*alpha);
				}
				else{
					weights(i,j) = FLT_MAX;
				}
			}
		}

		Eigen::MatrixXf P = PointCloud2Eigen(original);
		Eigen::MatrixXf Q = PointCloud2Eigen(shifted);
		Eigen::Matrix3f PQt;
		Eigen::Vector3f Pstar, Qstar;
		const Eigen::MatrixXf weightSum = weights.rowwise().sum();
		Eigen::ArrayXf temp = P;

		for(int n=0;n<cloudSize;n++){
			//PQt.setZero();
			Pstar = (weights.row(n)*P.transpose()).transpose()/weightSum(n);
			Qstar = (weights.row(n)*Q.transpose()).transpose()/weightSum(n);
			PQt = ((P.colwise()-Pstar).array().rowwise()*weights.row(n).array())
				     .matrix()*(Q.colwise()-Qstar).transpose();

			Eigen::JacobiSVD<Eigen::MatrixXf> svd(PQt, Eigen::ComputeFullU | Eigen::ComputeFullV);
			float scaling = svd.singularValues().trace()/(P.transpose()*P).trace();
			final->push_back(E2P((svd.matrixV()*svd.matrixU().transpose()*
			                  ((cloud->points[n]).getVector3fMap()-Pstar)*scaling*2+Qstar)));
			}
		return final;
		}
};

boost::shared_ptr<pcl::visualization::PCLVisualizer> simpleVis
  (typename pcl::PointCloud<pcl::PointXYZ>::ConstPtr cloud){
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer
	  (new pcl::visualization::PCLVisualizer ("3D Viewer"));
	viewer->setBackgroundColor (0, 0, 0);
	viewer->addPointCloud<pcl::PointXYZ> (cloud, "sample cloud");
	viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "sample cloud");
	//viewer->addCoordinateSystem (0.1);
	viewer->initCameraParameters ();
	return (viewer);
}

boost::tuple< pcl::PointCloud<pcl::PointXYZ>::Ptr,pcl::PointCloud<pcl::PointXYZ>::Ptr,
	pcl::PointCloud<pcl::PointXYZ>::Ptr >
	createPlateCloud(float Ra,int rc=100,int ac=100,float X = 0.0,float Y = 0.0,float Z = 0.0){
	pcl::PointCloud<pcl::PointXYZ>::Ptr out(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::PointCloud<pcl::PointXYZ>::Ptr center(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::PointCloud<pcl::PointXYZ>::Ptr rim(new pcl::PointCloud<pcl::PointXYZ>());
	Eigen::VectorXf T(ac,1);
	T.setLinSpaced(100,0,6.2830);
	Eigen::VectorXf R(rc,1);
	R.setLinSpaced(rc,0,Ra);
	for(int t=0;t<T.rows();t++){
		for(int r=0;r<R.rows();r++){
			if((r>=0 && r<=3)){
				center->push_back(pcl::PointXYZ(X+R(r,0)*sin(T(t,0)),Y+R(r,0)*cos(T(t,0)),Z));
			}
			else if(r>R.rows()-7){
				rim->push_back(pcl::PointXYZ(X+R(r,0)*sin(T(t,0)),Y+R(r,0)*cos(T(t,0)),Z));
			}
		 out->push_back(pcl::PointXYZ(X+R(r,0)*sin(T(t,0)),Y+R(r,0)*cos(T(t,0)),Z));
		}
	}
	return boost::tuple< pcl::PointCloud<pcl::PointXYZ>::Ptr,pcl::PointCloud<pcl::PointXYZ>::Ptr,
			pcl::PointCloud<pcl::PointXYZ>::Ptr >(out,center,rim);
}

int main(int argv,char *args[]){
	if(argv<2){
		std::cerr<<"Usage mls2 alpha"<<endl;
		return 1;
	}
	const float alpha = atof(args[1]);

	boost::tuple<pcl::PointCloud<pcl::PointXYZ>::Ptr,pcl::PointCloud<pcl::PointXYZ>::Ptr,
							  pcl::PointCloud<pcl::PointXYZ>::Ptr> temp =
								createPlateCloud(5,1000,1000);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = boost::tuples::get<0>(temp);
	pcl::PointCloud<pcl::PointXYZ>::Ptr center = boost::tuples::get<1>(temp);
	pcl::PointCloud<pcl::PointXYZ>::Ptr rim = boost::tuples::get<2>(temp);
	pcl::PointCloud<pcl::PointXYZ>::Ptr original(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr shifted(new pcl::PointCloud<pcl::PointXYZ>);
	const long cloudSize = cloud->size();

	original->insert(original->end(),rim->begin(),rim->end());
	original->insert(original->end(),center->begin(),center->end());

	for(pcl::PointCloud<pcl::PointXYZ>::iterator i = center->begin();i < center->end();i++){
		i->z-=2;
	}

	shifted->insert(shifted->end(),rim->begin(),rim->end());
	shifted->insert(shifted->end(),center->begin(),center->end());

  std::cout<<"Beginning deformation"<<std::endl;
	MLSDeform<pcl::PointXYZ> deformer(cloud,original,shifted,alpha);
	pcl::PointCloud<pcl::PointXYZ>::Ptr final(deformer.deform());
	std::cout<<"Done!"<<std::endl;

	boost::shared_ptr<pcl::visualization::PCLVisualizer> v = simpleVis(final);
	v->addPointCloud<pcl::PointXYZ>(original,"original");
	v->addPointCloud<pcl::PointXYZ>(shifted,"modified");
	while(!v->wasStopped()){
		v->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(10000));
	}
}
