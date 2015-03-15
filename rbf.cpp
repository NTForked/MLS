#include<iostream>
#include<vector>

#include<Eigen/Dense>
#include<pcl/common/common_headers.h>
#include<pcl/visualization/pcl_visualizer.h>

#include<boost/thread.hpp>
#include<boost/tuple/tuple.hpp>

boost::shared_ptr<pcl::visualization::PCLVisualizer> simpleVis(){
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer
	  (new pcl::visualization::PCLVisualizer ("3D Viewer"));
	viewer->setBackgroundColor (0, 0, 0);
	//viewer->addCoordinateSystem (0.1);
	viewer->initCameraParameters ();
	return (viewer);
}

float eucdist(const pcl::PointXYZ &a,const pcl::PointXYZ &b){
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
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

Eigen::MatrixXf RBFdeform(pcl::PointCloud<pcl::PointXYZ>::Ptr initial,pcl::PointCloud<pcl::PointXYZ>::Ptr final){
	int n = initial->size();
	Eigen::MatrixXf R(n,n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			R(i,j) = pow(eucdist(initial->points[i],initial->points[j]),3);
		}
	}
	Eigen::MatrixXf B(n,3);
	for(int i=0;i<n;i++){
		B(i,0) = final->points[i].x - initial->points[i].x;
		B(i,1) = final->points[i].y - initial->points[i].y;
		B(i,2) = final->points[i].z - initial->points[i].z;
	}
	return R.lu().solve(B);
}

int main(){
	boost::tuple<pcl::PointCloud<pcl::PointXYZ>::Ptr,pcl::PointCloud<pcl::PointXYZ>::Ptr,
							   pcl::PointCloud<pcl::PointXYZ>::Ptr> temp =
								createPlateCloud(5,1000,1000);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = boost::tuples::get<0>(temp);
	pcl::PointCloud<pcl::PointXYZ>::Ptr center = boost::tuples::get<1>(temp);
	pcl::PointCloud<pcl::PointXYZ>::Ptr rim = boost::tuples::get<2>(temp);
	pcl::PointCloud<pcl::PointXYZ>::Ptr original(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr shifted(new pcl::PointCloud<pcl::PointXYZ>);

	original->insert(original->end(),rim->begin(),rim->end());
	original->insert(original->end(),center->begin(),center->end());

	for(pcl::PointCloud<pcl::PointXYZ>::iterator i = center->begin();i < center->end();++i){
		i->z-=2;
	}

	shifted->insert(shifted->end(),rim->begin(),rim->end());
	shifted->insert(shifted->end(),center->begin(),center->end());
	std::cout<<"Begin"<<std::endl;
	RBFdeform(original,shifted);
	std::cout<<"Finish"<<std::endl;
}

