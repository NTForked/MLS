#include<iostream>
#include<vector>

#include<Eigen/Dense>
#include<pcl/common/common_headers.h>
#include<pcl/common/pca.h>
#include<pcl/common/centroid.h>
#include<pcl/common/distances.h>
#include<pcl/visualization/pcl_visualizer.h>

#include<boost/thread.hpp>
#include<boost/tuple/tuple.hpp>
#include<boost/lexical_cast.hpp>

boost::shared_ptr<pcl::visualization::PCLVisualizer> simpleVis(){
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer
	  (new pcl::visualization::PCLVisualizer ("3D Viewer"));
	viewer->setBackgroundColor (0, 0, 0);
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
		 out->push_back(pcl::PointXYZ(X+R(r,0)*sin(T(t,0)),Y+R(r,0)*cos(T(t,0)),Z+0.1*pow(R(r,0),3)));
		}
	}
	return boost::tuple< pcl::PointCloud<pcl::PointXYZ>::Ptr,pcl::PointCloud<pcl::PointXYZ>::Ptr,
			pcl::PointCloud<pcl::PointXYZ>::Ptr >(out,center,rim);
}

pcl::PointCloud<pcl::PointXYZ>::Ptr createLineCloud(){
	pcl::PointCloud<pcl::PointXYZ>::Ptr line(new pcl::PointCloud<pcl::PointXYZ>());
	for(int i=-100;i<100;i++){
		line->push_back(pcl::PointXYZ(i,0,0));
	}
	return line;
}

class HierachicalClustering{
	/*
	Use a top down approach to subdivide a point cloud.
	Recursively split a point cloud into clusters until:
		- The VariationalCriterion returns an number smaller than threshold
	*/
	public:
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> final;	
	class Plane{
	public:
		const Eigen::Matrix<float,3,1> normal;
		const Eigen::Matrix<float,3,1> centroid;
		Plane(Eigen::Matrix<float,3,1> normal,Eigen::Matrix<float,3,1> centroid):
													normal(normal),
													centroid(centroid){}

		bool query(pcl::PointXYZ &in){
			if( normal.transpose()*in.getVector3fMap() >= normal.transpose()*centroid){
				return true;
			}
			else
				return false;
		}
	};

	HierachicalClustering(pcl::PointCloud<pcl::PointXYZ>::Ptr &in){
		final.reserve(100);
		divide(in,final);
	}

	void divide(pcl::PointCloud<pcl::PointXYZ>::Ptr in,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> &final){
		if(in->size()<3)
			return;
		pcl::PCA<pcl::PointXYZ> pca;
		pca.setInputCloud(in);
		Eigen::MatrixXf axes = pca.getEigenValues();
		double planarity = axes(0)*axes(1)/axes(2)/axes.norm();
		//Do planarity analysis
		//std::cout<<planarity<<"  "<<in->size()<<std::endl;
		if(planarity<100000 && in->size()>100){
			Eigen::Matrix3f dirs = pca.getEigenVectors();
			Eigen::Matrix<float,4,1> _centroid;
			pcl::compute3DCentroid(*in,_centroid);
			Eigen::Matrix<float,3,1> centroid;
			for(int i=0;i<=2;i++){
				centroid(i) = _centroid(i)/_centroid(3); 
			}
			Eigen::Matrix<float,3,1> normal(dirs.col(0));
			// std::cout<<"normal="<<normal.transpose()<<std::endl;
			// std::cout<<"centroid="<<centroid.transpose()<<std::endl;
			// getchar();
			Plane sep(normal,centroid);
			pcl::PointCloud<pcl::PointXYZ>::Ptr left(new pcl::PointCloud<pcl::PointXYZ>());
			pcl::PointCloud<pcl::PointXYZ>::Ptr right(new pcl::PointCloud<pcl::PointXYZ>());
			for(pcl::PointCloud<pcl::PointXYZ>::iterator i = in->begin();i<in->end();++i){
				if(sep.query(*i)){
					left->push_back(*i);
				}
				else{
					right->push_back(*i);
				}
			}
		//std::cout<<"Calling left\nLength="<<left->size()<<std::endl;
		divide(left,final);
		//std::cout<<"Calling right\nLength="<<right->size()<<std::endl;
		divide(right,final);
		}
		else{
			final.push_back(in);
		}
	}

};

int main(){
	boost::tuple<pcl::PointCloud<pcl::PointXYZ>::Ptr,pcl::PointCloud<pcl::PointXYZ>::Ptr,
							   pcl::PointCloud<pcl::PointXYZ>::Ptr> temp =
								createPlateCloud(5,100,100);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = boost::tuples::get<0>(temp);
	HierachicalClustering a(cloud);
	boost::shared_ptr<pcl::visualization::PCLVisualizer> v = simpleVis();
	std::cout<<a.final.size()<<std::endl;
	srand(time(NULL));
	for(int i=0;i<a.final.size();i++){
		pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color(a.final[i], 
			(int)(rand()%255)+1, (int)(rand()%255)+1, (int)(rand()%255)+1);
		//pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color(cloud,255,255,255);
		v->addPointCloud(a.final[i],single_color,boost::lexical_cast<std::string>(i));
		v->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, boost::lexical_cast<std::string>(i));
	}
	// v->addPointCloud<pcl::PointXYZ> (a.final[0], "sample cloud");
	// v->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "sample cloud");
	v->initCameraParameters();
	std::cout<<"Done"<<std::endl;
	while(!v->wasStopped()){
		v->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(10000));
	}
}