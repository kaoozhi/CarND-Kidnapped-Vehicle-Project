/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::find_if;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  

  // Create random number generator
  std::default_random_engine gen;
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  double weight =1.0;

  // Create a normal distribution for first measurements given noises
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  // Assign generated states according to normal distribution to each particle
  for (int i = 0; i < num_particles; ++i) {    
    Particle sample;
    sample.id = 1;
    sample.x = dist_x(gen);
    sample.y = dist_y(gen);
    sample.theta = dist_theta(gen); 
    sample.weight = weight;
    particles.push_back(sample);
    weights.push_back(weight);
  }
  
  is_initialized = true;
  return;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  // Get measurement noises
  std::default_random_engine gen;
  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

 
  for (int i=0; i<particles.size();++i)
  {    
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;

    double x_final, y_final, theta_final;    
    
    // Motion model if yaw_rate is close to 0
    if (fabs(yaw_rate) <0.00001){
      x_final = x + velocity*delta_t*cos(theta);
      y_final = y + velocity*delta_t*sin(theta);
      theta_final = theta;
    }

    // Nominal motion model
    else
    {
      x_final= x + velocity/yaw_rate*(sin(theta + yaw_rate*delta_t)-sin(theta));
      y_final = y + velocity/yaw_rate*(cos(theta) - cos(theta + yaw_rate*delta_t));
      theta_final = theta + yaw_rate * delta_t;
    }

    // Create a normal distribution for particle's predicted positions given noises
    normal_distribution<double> noisy_x(x_final, std_x);
    normal_distribution<double> noisy_y(y_final, std_y);
    normal_distribution<double> noisy_theta(theta_final, std_theta);

    particles[i].x = noisy_x(gen);
    particles[i].y = noisy_y(gen);
    particles[i].theta = noisy_theta(gen);
    particles[i].theta = fmod(particles[i].theta, 2*M_PI);

  }
  return;
}

void ParticleFilter::dataAssociation(LandmarkObs &observation, 
                                     const vector<Map::single_landmark_s> &predicted_landmarks) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  double distance_min = std::numeric_limits<double>::max();
  int closest_id;
    for (int i=0; i<predicted_landmarks.size(); ++i)
    {
        Map::single_landmark_s landmark = predicted_landmarks[i];
        double distance = dist(double(landmark.x_f), double(landmark.y_f), observation.x, observation.y);
        
        if (distance<distance_min)
        {
            distance_min = distance;
            closest_id = landmark.id_i;
        }
    }
    observation.id = closest_id;
    
  return;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  double sig_x = std_landmark[0];
  double sig_y = std_landmark[1];

  // vector<double> update_weights; 
  for (int i=0; i<particles.size(); ++i){
      //get particle's coordinates
      double x_p = particles[i].x;
      double y_p = particles[i].y;
      double theta_p = particles[i].theta;


      vector<double> sense_x;
      vector<double> sense_y;
      vector<int> associations;
      double weight =1.0;

      // Filter lankmarks to keep only plausible ones within sensor range
      vector<Map::single_landmark_s> predicted_landmarks;
      for (int l=0; l<map_landmarks.landmark_list.size(); ++l)
      {
          Map::single_landmark_s current_landmark = map_landmarks.landmark_list[l];
          if ((fabs((x_p - current_landmark.x_f)) <= sensor_range) && (fabs((y_p - current_landmark.y_f)) <= sensor_range)) predicted_landmarks.push_back(current_landmark);
      }

      // Transform observations in MAP's coordinate with respect to particle position
      for (int j=0; j<observations.size(); ++j)
      {
        double x_c = observations[j].x;
        double y_c = observations[j].y;
        double x_m = x_p + cos(theta_p)*x_c - sin(theta_p)*y_c;
        double y_m = y_p + sin(theta_p)*x_c + cos(theta_p)*y_c;

        LandmarkObs transformed_observation;
        transformed_observation.id = observations[j].id;
        transformed_observation.x = x_m;
        transformed_observation.y = y_m;

        // Find nearest landmark in map and associate with observation
        dataAssociation(transformed_observation, predicted_landmarks);

        // Get associated map landmark's coordinates
        int id = transformed_observation.id;
        auto it = find_if(predicted_landmarks.begin(), predicted_landmarks.end(), [=](const Map::single_landmark_s & list) {return list.id_i == id;});

        // Caculate muiltivariable gaussian probability of transformed observations compared to associated map's landmark
        double mu_x = double(it->x_f);
        double mu_y = double(it->y_f);

        // Combine weights for current particle
        weight *= multiv_prob(sig_x, sig_y, x_m, y_m, mu_x, mu_y);

        // Append associations to currrent particle
        associations.push_back(transformed_observation.id);
        sense_x.push_back(transformed_observation.x);
        sense_y.push_back(transformed_observation.y);
      }

      // Update
      particles[i].weight = weight;
      weights[i]=weight;
      SetAssociations(particles[i],associations, sense_x, sense_y);
  }
  return;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

    std::default_random_engine gen;
    vector<Particle> resampled_particles;

    // Resampling Wheel method
    // int idx_max = particles.size()-1;
    // std::discrete_distribution<int> rand_idx(0, idx_max);
    // double w_max = *std::max_element(weights.begin(), weights.end());
    // std::uniform_real_distribution<double> rand_u(0, 2*w_max);
    // double beta = 0;
    // int idx = rand_idx(gen);
    // for (int i=0; i<particles.size(); ++i)
    // {
    //   double u = rand_u(gen);
    //   beta += u;
    //   while (weights[idx]< beta)
    //   {
    //       beta -= weights[idx];
    //       idx = fmod(idx+1, idx_max+1);
    //   }

    //   resampled_particles.push_back(particles[idx]);
    // }

    // Weighted resampling method
    std::discrete_distribution<int> rand_idx(weights.begin(), weights.end());
    for (int i=0; i<particles.size(); ++i)
    {
      resampled_particles.push_back(particles[rand_idx(gen)]);
    }
    particles = resampled_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}