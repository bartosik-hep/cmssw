//$Id: SprTrainedBagger.cc,v 1.6 2007/07/11 19:52:13 narsky Exp $

#include "PhysicsTools/StatPatternRecognition/interface/SprExperiment.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprTrainedBagger.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprUtils.hh"
#include "PhysicsTools/StatPatternRecognition/interface/SprDefs.hh"

#include <stdio.h>
#include <cassert>

using namespace std;


SprTrainedBagger::SprTrainedBagger(const std::vector<
		   std::pair<const SprAbsTrainedClassifier*,bool> >& 
		   trained, bool discrete) 
  : 
  SprAbsTrainedClassifier(),
  trained_(trained),
  discrete_(discrete)
{
  assert( !trained_.empty() );
  this->setCut(SprUtils::lowerBound(0.5));
}


SprTrainedBagger::SprTrainedBagger(const SprTrainedBagger& other)
  :
  SprAbsTrainedClassifier(other),
  trained_(),
  discrete_(other.discrete_)
{
  for( int i=0;i<other.trained_.size();i++ )
    trained_.push_back(pair<const SprAbsTrainedClassifier*,bool>
		       (other.trained_[i].first->clone(),true));
}


double SprTrainedBagger::response(const std::vector<double>& v) const
{
  // init
  double r = 0;

  // discrete/continuous
  if( discrete_ ) {
    int out = 0;
    for( int i=0;i<trained_.size();i++ )
      out += ( trained_[i].first->accept(v) ? 1 : -1 );
    r = out;
    r /= 2.*trained_.size();
    r += 0.5;
  }
  else {
    for( int i=0;i<trained_.size();i++ )
      r += trained_[i].first->response(v);
    r /= trained_.size();
  }

  // exit
  return r;
}


void SprTrainedBagger::destroy()
{
  for( int i=0;i<trained_.size();i++ ) {
    if( trained_[i].second )
      delete trained_[i].first;
  }
}


void SprTrainedBagger::print(std::ostream& os) const
{
  os << "Trained Bagger " << SprVersion << endl;
  os << "Classifiers: " << trained_.size() << endl;
  for( int i=0;i<trained_.size();i++ ) {
    os << "Classifier " << i 
       << " " << trained_[i].first->name().c_str() << endl;
    trained_[i].first->print(os);
  }
}


bool SprTrainedBagger::generateCode(std::ostream& os) const 
{ 
  // generate weak classifiers
  for( int i=0;i<trained_.size();i++ ) { 
    string name = trained_[i].first->name();
    os << " // Classifier " << i  
       << " \"" << name.c_str() << "\"" << endl; 
    if( !trained_[i].first->generateCode(os) ) {
      cerr << "Unable to generate code for classifier " << name.c_str() 
	   << endl;
      return false;
    }
    if( i < trained_.size()-1 ) os << endl; 
  }

  // exit
  return true; 
} 

