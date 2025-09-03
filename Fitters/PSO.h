#pragma once

// C++ includes
#include <algorithm>
#include <functional>
#include <numeric>

// MaCh3 includes
#include "Fitters/LikelihoodFit.h"


/// @brief Class particle - stores the position, velocity and personal best
/// With functions which move particle and update velocity
/// @author Emily Ip
/// @author Mark Scott
/// @date 24/2/2023
class particle{
  public:
      /// @brief Constructor
      particle(const std::vector<double>& pos, const std::vector<double>& vel) : position(pos), velocity(vel){};
      /// @brief Destructor
      virtual ~particle() {};

      void set_position(const std::vector<double>& new_position) {
        position = new_position;
      };

      std::vector<double> get_position() {
        return position;
      };

      void set_personal_best_position(const std::vector<double>& new_pos){
        personal_best_position = new_pos;
      };

      std::vector<double> get_personal_best_position(){
        return personal_best_position;
      };

      void set_personal_best_value(const double new_val){
        personal_best_value = new_val;
      };

      double get_personal_best_value(){
        return personal_best_value;
      };

      std::vector<double> get_velocity(){
        return velocity;
      };

      void set_velocity(const std::vector<double>& new_velocity){
        velocity = new_velocity;
      };

      double get_value(){
        return curr_value;
      };
      void set_value(const double val){
        curr_value = val;
      };

  private:
    std::vector<double> position;
    std::vector<double> velocity;
    double personal_best_value;
    double curr_value;
    std::vector<double> personal_best_position;
};

 /// @brief Class PSO, consist of a vector of object Class Particle and global best
 /// Takes in the size (number of particle) and number of iteration
 /// functions includes: finding global best, updating velocity, actual minimisation function
/// @author Emily Ip
/// @author Mark Scott
/// @date 24/2/2023
class PSO : public LikelihoodFit {
  public:
    /// @brief constructor
    PSO(manager * const fitMan);
    /// @brief Destructor
    virtual ~PSO() {};

    particle* get_best_particle(){
        return best_particle;
    }
    void set_best_particle(particle* n){
      best_particle = n;
    }

    std::vector<std::vector<double>> bisection(const std::vector<double>& position, const double minimum,
                                               const double range, const double precision);
    std::vector<std::vector<double>> calc_uncertainty(const std::vector<double>& position, const double minimum);
    void init();
    void uncertainty_check(const std::vector<double>& previous_pos);
    void run();
    void WriteOutput();
    /// @brief Actual implementation of PSO Fit algorithm
    void RunMCMC() override;
    double CalcChi2(const double* x) override;
    /// @brief Evaluates the Rastrigin function for a given parameter values.
    double rastriginFunc(const double* x);
    double swarmIterate();

    std::vector<double> vector_multiply(std::vector<double> velocity, const double mul){
      // std::bind1st deprecated since C++11, removed in c++17
      // transform(velocity.begin(),velocity.end(),velocity.begin(),std::bind1st(std::multiplies<double>(),mul));
      std::transform(velocity.begin(), velocity.end(), velocity.begin(),
                      std::bind(std::multiplies<double>(), mul, std::placeholders::_1));
      return velocity;
    };

    std::vector<double> vector_add(const std::vector<double>& v1, const std::vector<double>& v2){
      std::vector<double> v3;
      transform(v1.begin(), v1.end(), v2.begin(), back_inserter(v3), std::plus<double>());
      return v3;
    };
    std::vector<double> vector_subtract(const std::vector<double>& v1, const std::vector<double>& v2){
      std::vector<double> v3 ;
      transform(v1.begin(), v1.end(), v2.begin(), back_inserter(v3), std::minus<double>());
      return v3;
    };
    std::vector<double> three_vector_addition(std::vector<double> vec1,
                                              const std::vector<double>& vec2,
                                              const std::vector<double>& vec3) {
      for (size_t i = 0; i < vec1.size(); ++i) {
        vec1[i] += vec2[i] + vec3[i];
      }
      return vec1;
    };
    std::vector<double> four_vector_addition(std::vector<double> vec1, const std::vector<double>& vec2,
                                             const std::vector<double>& vec3, const std::vector<double>& vec4){
      for (size_t i = 0; i < vec1.size(); ++i) {
        vec1[i] += vec2[i] + vec3[i] + vec4[i];
      }
      return vec1;
    };

    double CalcChi(std::vector<double> x){
      double* a = &x[0];
      return CalcChi2(a);
    };

  inline std::string GetName() const override {return "PSO";};

  private:
    particle* best_particle;
    double fBestValue;
    std::vector<double> prior;
    std::vector<bool> fixed;
    std::vector<double> ranges_max;
    std::vector<double> ranges_min;
    std::vector<particle*> system;
    double fInertia;
    double fOne;
    double fTwo;
    double fConvergence;
    int fIterations;
    double fConstriction;
    std::vector<std::vector<double> > uncertainties;

    int fParticles;
    std::vector<double*> paramlist;
    constexpr static const int kMaxParticles = 10000;
    double vel[kMaxParticles];
    double* par;

    int fDim;
};

