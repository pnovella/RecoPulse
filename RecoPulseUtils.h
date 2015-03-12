

#ifndef RecoPulseUtils_hh
#define RecoPulseUtils_hh

#include<string>
#include<vector>
#include<cmath>
#include <iostream>


static const double NOPED = -1;

static const int NOTIME = -1;

static const int NOCAL = 0;



// class messenger{


// private:
  
//   size_t level_;
  
//   std::string name_;

// public:
  
  
//   messenger(){name_="RPMessenger",level_ = 0;};

//   messenger(size_t clevel){name_="RPMessenger",level_ = clevel;}
  
//   messenger(std::string name, size_t clevel=0){name_=name;level_ = clevel;}

//   void set_level(size_t clevel){level_ = clevel;}

//   void set_name(std::string name){name_ = name;}
 
//   size_t level() const {return level_;} 


//   inline void message (const std::string message, 
// 		       size_t clevel,std::string h="") const
//   {
//     if (!h.empty()) h += ": ";
 
//     if (clevel <= level_)
//       cout<<"+++ "<<name_<<"("<<level_<<":"<<clevel<<"): "<<h
// 	   << message << endl;
//   }

//   template <class T> inline
//   void message(const std::string message, const T& d, 
// 	       size_t clevel,std::string h="") const
//   {
//     if (!h.empty()) h += ": ";

//     if (clevel <= level_)
//       cout<<"+++ "<<name_<<"("<<level_<<":"<<clevel<<"): "<< h<<
// 	message << " " << d << endl;
//   }


//   template <class A, class B> inline
//   void message(const std::string message, const A& d1, 
// 	       const B& d2, size_t clevel,std::string h="") const
//   {
//     if (!h.empty()) h += ": ";

//     if (clevel <= level_)
//       cout<<"+++ "<<name_<<"("<<level_<<":"<<clevel<<"): "<<h
// 	   << message << " " << d1 << " " << d2 << endl;
//   }

//   inline void warning (const std::string message, 
// 		       size_t clevel) const
//   {
//     this->message(message,clevel,"WARNING");
//   }

//   template <class T> inline
//   void warning(const std::string message, const T& d, 
// 	       size_t clevel) const
//   {
//     this->message(message,d,clevel,"WARNING");
//   }


//   template <class A, class B> inline
//   void warning(const std::string message, const A& d1, 
// 	       const B& d2, size_t clevel) const
//   {
//     this->message(message,d1,d2,clevel,"WARNING");
//   }
  
//    inline void error (const std::string message, 
// 		       size_t clevel) const
//   {
//     this->message(message,clevel,"ERROR");
//   }

//   template <class T> inline
//   void error(const std::string message, const T& d, 
// 	       size_t clevel) const
//   {
//     this->message(message,d,clevel,"ERROR");
//   }


//   template <class A, class B> inline
//   void error(const std::string message, const A& d1, 
// 	       const B& d2, size_t clevel) const
//   {
//     this->message(message,d1,d2,clevel,"ERROR");
//   }

// };


class spline3{
  
private:
  
  std::vector<double> _nodes; std::vector<double> _values;

  std::vector<double> _acoefs, _bcoefs, _ccoefs,_dcoefs;
  
  bool _status;

private:

  double integral(double,size_t);

  int getIndex(double x);

public:
  
  spline3();
  
  spline3(const std::vector<unsigned short>&);

  spline3(const std::vector<double>&);

  spline3(const std::vector<double>&,
	  const std::vector<double>&);
  
  virtual ~spline3(){}
  
  void gen();

  void setNodes(const std::vector<double>&);

  void setValues(const std::vector<double>&);
  

public:

  
  double eval(double x);

  
  double values(size_t i) const {return _values[i];}

  double nodes(size_t i) const {return _nodes[i];}
  
  size_t size() const {return _nodes.size();}
  
  bool status() const {return _status;}
  
  double integrate(size_t min=0,size_t max=0);

};

class peakGatherer{

private:
  
  std::vector<unsigned short> _prof;

  std::vector<int> _T0s, _T1s;
  
  bool _status;
  

public:
  
  peakGatherer(){_status=false;}

  peakGatherer(const std::vector<unsigned short>& prof)
  {_prof=prof; _status=false;}


  virtual ~peakGatherer(){};
  
  //! get peaks with positive charge
  void getPeaks(double thr,size_t nsamp=3);
  
  //! get bi-polar peaks 
  void getBPPeaks(double thr,size_t nsamp=3);
  
  std::vector<int> getT0s() const {return _T0s;}

  std::vector<int> getT1s() const {return _T1s;}

  bool status() const {return _status;}

};


class RPPulse {

private:

  std::vector<unsigned short> _prof;

  std::vector<unsigned short> _subProf;

  // interpolation
  spline3 _spline;

  // peak search
  peakGatherer _peaks;
  
  // integrated bins
  std::vector<bool> _integrated;

public:
  
  RPPulse(const std::vector<unsigned short>& prof){
    
    _prof = prof; 
    
    //_spline = spline3(_prof);

    _peaks = peakGatherer(_prof);
    
    _integrated = std::vector<bool>(_prof.size(),false);

  }
  
  RPPulse(const RPPulse& pulse){
    
    _prof =std::vector<unsigned short>(pulse.getProfile());

    //_spline = spline3(_prof);

    _peaks = peakGatherer(_prof);

    _integrated = std::vector<bool>(_prof.size(),false);

  }
  
  virtual ~RPPulse(){}
  
  void genSpline(){_spline = spline3(_prof); _spline.gen();}
  
  void getPeaks(double thr,size_t nsamp=10){_peaks.getPeaks(thr,nsamp);}

  void getBPPeaks(double thr,size_t nsamp=10){_peaks.getBPPeaks(thr,nsamp);}
  
  void filter(double ped, double thr,size_t nsamples=3);

  //test
  bool hasPeak(double ped, double thr, size_t nsamples,
		    int t0, int t1);
  
  size_t getContSamples(double ped,size_t t0=0, size_t t1=0);

  size_t getPeakWidth(double ped, double thr,int t0=0, int t1=0);

public:
  
  spline3 getSpline(){return _spline;}
  
  peakGatherer getPeakGatherer(){return _peaks;}

  const std::vector<unsigned short>& getProfile() const {return _prof;} 

  size_t size() const {return _prof.size();}
  
  const std::vector<unsigned short>& getSubProfile(double ped);
  
  void smoothPed(double,size_t,size_t);

  void integrated(size_t,size_t);

  bool isIntegrated(size_t i) const {return _integrated[i];}

};




#endif
/** @} */ // end of doxygen group RecoPulse
