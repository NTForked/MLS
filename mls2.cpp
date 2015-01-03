#include<iostream>
#include<vector>

#include<cmath>
#include<cstdlib>
#include<cfloat>

#include<Eigen/Dense>

#include<pcl/common/common_headers.h>
#include<pcl/common/centroid.h>
#include<pcl/visualization/pcl_visualizer.h>

#include<boost/thread/thread.hpp>
#include<boost/tuple/tuple.hpp>

typedef Eigen::MatrixXf Matrix;
typedef Eigen::VectorXf Vector;

typedef pcl::PointXYZ PointPCL;
typedef pcl::PointCloud<PointPCL> PointCloud;
typedef PointCloud::Ptr Ptr;
typedef pcl::visualization::PCLVisualizer Visualizer;

using std::cout;
using std::endl;
using boost::tuple;
using boost::tuples::get;

template<class PointT>
class MovingLeastSquaresDeformation{
private:
	const pcl::PointCloud<PointT> cloud;
	const pcl::PointCloud<PointT> control;

	const float eucdist(const PointPCL &a,const PointPCL &b){
		return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
	}

	inline PointPCL E2P(const Matrix &in){
		return PointPCL(in(0,0),in(1,0),in(2,0));
	}

	inline Eigen::Matrix<float,3,1> P2E(const PointPCL &in){
		Eigen::Matrix<float,3,1> out;
		out<<in.x,in.y,in.z;
		return out;
	}

	Matrix PointCloud2Eigen(Ptr cloud){
		Matrix out(3,cloud->size());
		int j=0;
		for(PointCloud::iterator i=cloud->begin();i<cloud->end();i++){
			Eigen::MatrixXf temp(3,1);
			temp<<i->x,i->y,i->z;
			out.col(j++) = temp;
		}
		return out;
	}

	boost::shared_ptr<Visualizer> simpleVis (PointCloud::ConstPtr cloud){
		boost::shared_ptr<Visualizer> viewer (new Visualizer ("3D Viewer"));
		viewer->setBackgroundColor (0, 0, 0);
		viewer->addPointCloud<PointPCL> (cloud, "sample cloud");
		viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "sample cloud");
		//viewer->addCoordinateSystem (0.1);
		viewer->initCameraParameters ();
		return (viewer);
	}

	tuple<Ptr,Ptr,Ptr> createPlateCloud(float Ra,int rc=100,int ac=100,float X = 0.0,float Y = 0.0,float Z = 0.0){
		Ptr out(new PointCloud);
		Ptr center(new PointCloud);
		Ptr rim(new PointCloud);
		Vector T(ac,1);
		T.setLinSpaced(100,0,6.2830);
		Vector R(rc,1);
		R.setLinSpaced(rc,0,Ra);
		for(int t=0;t<T.rows();t++){
			for(int r=0;r<R.rows();r++){
				if((r>=0 && r<=3)){
					center->push_back(PointPCL(X+R(r,0)*sin(T(t,0)),Y+R(r,0)*cos(T(t,0)),Z));
				}
				else if(r>R.rows()-7){
					rim->push_back(PointPCL(X+R(r,0)*sin(T(t,0)),Y+R(r,0)*cos(T(t,0)),Z));
				}
				else{
					out->push_back(PointPCL(X+R(r,0)*sin(T(t,0)),Y+R(r,0)*cos(T(t,0)),Z));
				}
			}
		}
		return tuple<Ptr,Ptr,Ptr>(out,center,rim);
	}
	
public:
	MovingLeastSquaresDeformation(const pcl::PointCloud<PointT> &cloud_,const pcl::PointCloud<PointT> &control_){

	}
	void deform(){

	}
};

int main(int argv,char *args[]){
	if(argv<2){
		std::cerr<<"Usage mls2 alpha"<<endl;
		return 1;
	}
	const float alpha = atof(args[1]);

		std::cerr<<"Creating cloud"<<endl;
	tuple<Ptr,Ptr,Ptr> temp = createPlateCloud(5);
	const Ptr cloud = get<0>(temp);
	const Ptr center = get<1>(temp);
	const Ptr rim = get<2>(temp);
	Ptr original(new PointCloud);
	Ptr shifted(new PointCloud);
	Ptr final(new PointCloud);
	const long cloudSize = cloud->size();
	final->reserve(cloudSize);

	original->insert(original->end(),rim->begin(),rim->end());
	original->insert(original->end(),center->begin(),center->end());

	for(PointCloud::iterator i = center->begin();i < center->end();i++){
		i->z-=2;
		}

	shifted->insert(shifted->end(),rim->begin(),rim->end());
	shifted->insert(shifted->end(),center->begin(),center->end());

	Matrix w(cloud->size(),original->size());
	bool flags[cloud->size()];

	for(int i=0;i<cloud->size();i++){
		flags[i] = true;
		for(int j=0;j<original->size();j++){
			if(eucdist(cloud->points[i],original->points[j])>0.0001){
				w(i,j)=pow(eucdist(cloud->points[i],original->points[j]),-2*alpha);
			}
			else{
				w(i,j)=FLT_MAX;
				flags[i]=false;
			}
		}
	}

		std:cerr<<"Finding optimal transformations.."<<endl;
	Matrix P = PointCloud2Eigen(original);
	Matrix Q = PointCloud2Eigen(shifted);
	Eigen::Matrix3f PQt;

	for(int n=0;n<cloud->size();n++){
		PQt.setZero();
		if(true){
			for(int i=0;i<original->size();i++){
				PQt += w(n,i)*P.col(i)*(Q.col(i).transpose());
			}

			Eigen::Vector3f Pstar, Qstar;
			const float weightSum = w.row(n).sum();
			for(int i=0;i<original->size();i++){
				Pstar+=w(n,i)*P.col(i);
				Qstar+=w(n,i)*Q.col(i);
			}
			Pstar/=weightSum;
			Qstar/=weightSum;

			Eigen::JacobiSVD<Matrix> svd(PQt, Eigen::ComputeFullU | Eigen::ComputeFullV);
			float scaling = svd.singularValues().trace()/(P.transpose()*P).trace();
			final->push_back(E2P((svd.matrixV()*svd.matrixU().transpose()*(P2E(cloud->points[n])-Pstar)*scaling/(1-weightSum/100)+Qstar)));
		}
	}

	boost::shared_ptr<Visualizer> v = simpleVis(final);
	v->addPointCloud<PointPCL>(original,"original");
	v->addPointCloud<PointPCL>(shifted,"modified");
	while(!v->wasStopped()){
		v->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(10000));
		}
	}
