/*
 *  bcpnn_connection.h
 *
 *  Written by Philip Tully
 *
 */

#ifndef BCPNN_CONNECTION_H
#define BCPNN_CONNECTION_H

/* BeginDocumentation
  Name: bcpnn_synapse - Synapse type for incremental, Bayesian spike-timing 
   dependent plasticity.

  Description:
   bcpnn_synapse is a connector to create synapses with incremental, Bayesian 
   spike timing dependent plasticity.

   tau_i	double - Primary trace presynaptic time constant
   tau_j	double - Primary trace postsynaptic time constant
   tau_e	double - Secondary trace time constant
   tau_p	double - Tertiarty trace time constant
   p_i		double - \
   p_j		double -  >- these 3 initial conditions determine weight, i.e. log(p_ij/(p_i * p_j)).
   p_ij		double - /
   K_		double - Print-now signal // Neuromodulation. Turn off learning, K = 0.
   fmax_        double - Frequency assumed as maximum firing, for match with abstract rule
   epsilon_     double - lowest possible probability of spiking, e.g. lowest assumed firing rate
   bias_        double - ANN interpretation. Only calculated here to demonstrate match to rule. 
                         Will be eliminated in future versions, where bias will be calculated postsynaptically
   gain_    double - Coefficient to scale weight as conductance, can be zero-ed out

  Transmits: SpikeEvent
   
  References:
   [1] Wahlgren and Lansner (2001) Biological Evaluation of a Hebbian-Bayesian
       learning rule. Neurocomputing, 38-40, 433-438

   [2] Bergel, Transforming the BCPNN Learning Rule for Spiking Units to a
       Learning Rule for Non-Spiking Units (2010). KTH Masters Thesis.

  FirstVersion: November 2011
  CurrentVersion: March 2012
  Author: Philip Tully
          tully@csc.kth.se
  SeeAlso: synapsedict, stdp_synapse, tsodyks_synapse, static_synapse
*/

/* for Debugging */
#include <iostream>
using namespace std;

#include "connection_het_wd.h"
#include "archiving_node.h"
#include "generic_connector.h"
#include <cmath>

namespace mynest
{
  class BCPNNConnection : public nest::ConnectionHetWD
  {
    public:
      /* Default Constructor. Sets default values for all parameters. Needed by GenericConnectorModel. */
      BCPNNConnection();

      /* Copy constructor. Needs to be defined properly in order for GenericConnector to work. */
      BCPNNConnection(const BCPNNConnection &);

      /* Default Destructor. */
      ~BCPNNConnection() {}

      void check_connection(nest::Node & s, nest::Node & r, nest::port receptor_type, nest::double_t t_lastspike);

      /* Get all properties of this connection and put them into a dictionary. */
      void get_status(DictionaryDatum & d) const;

      /* Set properties of this connection from the values given in dictionary. */
      void set_status(const DictionaryDatum & d, nest::ConnectorModel &cm);

      /* Set properties of this connection from position p in the properties array given in dictionary. */
      void set_status(const DictionaryDatum & d, nest::index p, nest::ConnectorModel &cm);

      /* Create new empty arrays for the properties of this connection in the given dictionary. It is assumed 
         that they do not exist before. */
      void initialize_property_arrays(DictionaryDatum & d) const;

      /* Append properties of this connection to the given dictionary. If the dictionary is empty, new arrays 
         are created first. */
      void append_properties(DictionaryDatum & d) const;

      /* Send an event to the receiver of this connection.  */
      void send(nest::Event& e, nest::double_t t_lastspike, const nest::CommonSynapseProperties &cp);

      /* Overloaded for all supported event types. */
      using nest::Connection::check_event;
      void check_event(nest::SpikeEvent&) {}

    private:
      /* data members of each connection */
      nest::double_t yi_;
      nest::double_t yj_;

      nest::double_t taui_;
      nest::double_t tauj_;
      nest::double_t taue_;
      nest::double_t taup_;

      nest::double_t epsilon_;
      nest::double_t K_;
      nest::double_t bias_;
      nest::double_t fmax_;
      nest::double_t gain_;

      nest::double_t zi_;
      nest::double_t zj_;
      nest::double_t ei_;
      nest::double_t ej_;
      nest::double_t eij_;
      nest::double_t pi_;
      nest::double_t pj_;
      nest::double_t pij_;
      nest::double_t t_k_;
	  std::vector<nest::double_t> times_k_changed;
	  std::vector<nest::double_t> post_spiketimes;
	  std::vector<nest::double_t> K_values_;
  }; /* of class BCPNNConnection */

  inline 
  void BCPNNConnection::check_connection(nest::Node & s, nest::Node & r, nest::port receptor_type, nest::double_t t_lastspike)
  {
    nest::ConnectionHetWD::check_connection(s, r, receptor_type, t_lastspike);

    // For a new synapse, t_lastspike contains the point in time of the last spike.
    // So we initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay]
    // which increases the access counter for these entries.
    // At registration, all entries' access counters of history[0, ..., t_last_spike - dendritic_delay] will be 
    // incremented by the following call to Archiving_Node::register_stdp_connection().
    // See bug #218 for details.
    r.register_stdp_connection(t_lastspike - nest::Time(nest::Time::step(delay_)).get_ms());
  }

  /* Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   * \param t_lastspike Time point of last spike emitted 
  
   note: every time this method is called by an outside function, a presynaptic
       event has occured and is being transmitted to the postsynaptic side. */

  inline
  void BCPNNConnection::send(nest::Event& e, nest::double_t t_lastspike, const nest::CommonSynapseProperties &)
  {
    nest::double_t t_spike = e.get_stamp().get_ms();  /* time stamp of current spike event */
    nest::double_t dendritic_delay = nest::Time(nest::Time::step(delay_)).get_ms();    /* delay from dendrite -> soma */
    nest::double_t resolution = nest::Time::get_resolution().get_ms();               /* nest.GetKernelStatus('resolution') simulation timestep */
    nest::int_t spike_width = 10;                     /* assume spike width of 1ms, resolution is 0.1 so mult by 10 */    
    nest::double_t spike_height = 1000.0 / fmax_;     /* normalizing to match this spiking rule to abstract = 1000/FMAX (Hz)*/
    nest::int_t counter = 0;                          /* ensuring traces reverberate for duration of the spike width */
    nest::double_t min_weight = epsilon_/std::pow(0.5 ,2);         /* theoretical minimum weight = epsilon/(0.5*0.5) */

//      cout << "t_lastspike: " << t_lastspike << " t_spike: " << t_spike << " K_values.size(): " << K_values_.size() << endl;

    /*STEP ONE: Get all timings of pre and postsynaptic spikes. Post store in dynamically allocated array */

    /* get spike history in relevant range (t1, t2] from post-synaptic neuron */
    std::deque<nest::histentry>::iterator start;
    std::deque<nest::histentry>::iterator finish;
  
    /* Initially read the history(t_last_spike - dendritic_delay, ...,  T_spike-dendritic_delay] which increases the 
       access counter for these entries. At registration, all entries' access counters of history[0, ..., 
       t_last_spike - dendritic_delay] have been incremented by Archiving_Node::register_stdp_connection(). 
       See bug #218 for details. */
    target_->get_history(t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish);

    /* For spike order pre-post, if dopamine present facilitate else depress.
       Pre  spikes: |       |  t_lastpike is the last pre spike and t_spike is the current pre spike
       Post spikes    | ||		 start is a pointer to the first post spike in the interval between the
       two pre spikes. It is then iterated until the last post spike in the interval */


    while (start != finish)  {/* loop until you get to last post spike */
		post_spiketimes.push_back(start->t_);
		start++;
    } /* of while */
    /* STEP TWO: Consider the presynaptic firing window, delta t resolution, and update the traces */
    
    /* nest stores with ms precision the timing of the spike. */
    /* the following loop iterates through the presynaptic spike difference window */
    counter = 0;
	nest::int_t number_iterations = (nest::int_t)((t_spike - t_lastspike)/resolution);
	// This if+else is to account for the case when K_ is initialized with 1 and K_ is changed via set_status /before/ the simulation starts
	std::vector<nest::double_t> K_vec (number_iterations, K_values_.front());
	if((nest::int_t)t_lastspike == 0) {
		K_vec = std::vector<nest::double_t> (number_iterations, K_values_.back());
	}

	std::vector<nest::double_t>::iterator post_it = post_spiketimes.begin(); 
	std::vector<nest::double_t>::iterator time_it = times_k_changed.end();
	std::vector<nest::double_t>::iterator K_it = K_values_.end();
	if (K_values_.size() > 1) { // only if a set_status has ever been called --> TODO
		if (times_k_changed.back() >= t_lastspike){ // if K has changed since the last pre-synaptic spike has occured
			// --> find out when changes occured and set the K_vec values right
			K_it--; // move iterator one element to the left, because .end() returns iterator just past the last element
			time_it--;
			nest::int_t idx_first = (nest::int_t) ((t_spike - t_lastspike) / resolution);
			nest::int_t idx_second;
			while (*time_it > t_lastspike){
				idx_second = (nest::int_t) ((*time_it - t_lastspike)/ resolution);
//                n_constant_K = (nest::int_t) ((idx_first - idx_second) / resolution);
//                cout << "New spike at *time_it: " << *time_it << " t_lastspike: " << t_lastspike << endl;
				for (nest::int_t i_k=idx_first-1; i_k >= idx_second; --i_k) { 
					// alternative implementation would be to create a temporary vector with the correct K-values and assign this to K_vec
//                    cout << "debug size K_vec = " << K_vec.size() << " i_k = " << i_k << " idx_first: " << idx_first << " idx_second: " << idx_second << " times_k_changed.size=" << times_k_changed.size() << " t_lastspike=" << t_lastspike << endl;
					K_vec.at(i_k) = *K_it;
				}
				idx_first = idx_second;
				time_it--;
				K_it--;
			} // end of while
//            K_ = *(K_values_.end()); // update the private K_ value
		}
	}
	K_values_.clear();
	K_values_.push_back(K_);
	times_k_changed.clear();
	times_k_changed.push_back(*time_it);


	nest::int_t idx_spike;

	/* Create a vector to represent the post spikes as a trace */
	std::vector<nest::double_t> post_active (number_iterations + spike_width, 0.);

//    std::vector<nest::double_t double> bar (5,0);
//    cout << "weight at send t_spike " << t_spike << " number_iterations: " << number_iterations << " t_lastspike: " << t_lastspike << endl;
	cout << "START PRINTOUT" << endl;
//    cout << t_spike << "\t" << number_iterations << "\t" << t_lastspike << endl;
    for (nest::int_t timestep = 0; timestep < number_iterations; timestep++)
    {
		/* CASE: Default. Neither Pre nor Post spike. */
		yi_ = 0.0; 
		yj_ = 0.0;

		/* CASE: Pre without (*OR WITH post) spike - synchronous events handled automatically. */
//        if(timestep < spike_width && (nest::int_t)t_lastspike != 0) {
		if(timestep == 0 && t_lastspike != 0.) {
			yi_ = spike_height * spike_width;
		}

		// if you have any post spike at all
		if (post_spiketimes.size() > 0) { 
			if (post_it != post_spiketimes.end()) { 
				if (timestep == (nest::int_t)((*post_it) - t_lastspike) / resolution){
					yj_ = spike_height * spike_width;
					post_it++;
				}
			}
		}

		/* Primary synaptic traces. Noise - commented out*/
		zi_ += (yi_ - zi_ + epsilon_ /*+ (0.01 + (double)rand() / RAND_MAX * (0.05 - 0.01))*/ ) * resolution / taui_;
//        zj_ += (spike_height * post_active.at(timestep) - zj_ + epsilon_ /*+ (0.01 + (double)rand() / RAND_MAX * (0.05 - 0.01))*/ ) * resolution / tauj_;
		zj_ += (yj_ - zj_ + epsilon_ /*+ (0.01 + (double)rand() / RAND_MAX * (0.05 - 0.01))*/ ) * resolution / tauj_;

		/* Secondary synaptic traces */
		ei_  += (zi_ - ei_) * resolution / taue_;
		ej_  += (zj_ - ej_) * resolution / taue_;
		eij_ += (zi_ * zj_ - eij_) * resolution / taue_;

		/* Tertiary synaptic traces. Commented is from Wahlgren paper. */
		//    pi_  += K_ * (ei_ - pi_) * resolution / taup_/* * eij_*/;
		//    pj_  += K_ * (ej_ - pj_) * resolution / taup_/* * eij_*/;
		//    pij_ += K_ * (eij_ - pij_) * resolution / taup_/* * eij_*/;
		pi_  += K_vec[timestep] * (ei_ - pi_) * resolution / taup_;
		pj_  += K_vec[timestep] * (ej_ - pj_) * resolution / taup_;
		pij_ += K_vec[timestep] * (eij_ - pij_) * resolution / taup_;

		weight_ = gain_ * std::log(pij_ / (pi_ * pj_)) /*- std::log(min_weight)*/;
		cout << timestep << "\t" << zi_ << "\t" << zj_ << "\t" << ei_ << "\t" << "\t" << ej_ << "\t" << eij_ << "\t" << pi_ << "\t" << pj_ << "\t" << pij_ << "\t" << weight_ << "\t" << t_spike << "\t" << yj_ << "\t" << yi_ << endl;
    } /* of for */
	cout << "END PRINTOUT" << endl;

    /* Update the weight & bias before event is sent. Use commented normalization to 
       implement soft weight bounds, this way the weight will never go below 0 because
       you push all weights up by the most negative weight possible. */
    bias_ = std::log(pj_);
    weight_ = gain_ * (std::log(pij_ / (pi_ * pj_)) /*- std::log(min_weight) */);

    /* STEP THREE. Implement hard weight bounds. NOTE if using above normalization, weights
                   are soft-bounded above zero already. */
    /*weight_ = (weight_ < 0) ? weight_ : 0.0;
      nest::double_t Wmax = ...;
      weight_ = (weight_ > Wmax) ? weight_ : Wmax;*/

    /* Send the spike to the target */
    e.set_receiver(*target_);
    e.set_weight(weight_);
    e.set_delay(delay_);
    e.set_rport(rport_);
    e();

    /* final clean up */
	// before clearing it, remember the last spike, for the case that the post synaptic spike is just before the next presynaptic spike
//    nest::double_t t_temp = post_spiketimes.back();
	post_spiketimes.clear();
//    post_spiketimes.push_back(0.);
//    cout << "after send, post_spiketimes holds: " << post_spiketimes.front() << endl;

//    K_vec.clear();

  } /* of BCPNNConnection::send */

} /* of namespace mynest */
#endif /* of #ifndef BCPNN_CONNECTION_H */

