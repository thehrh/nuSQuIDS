#ifndef __nuSQUIDNSI_H
#define __nuSQUIDNSI_H

#include "nuSQuIDS.h"

namespace nusquids{


	class nuSQUIDSNSI: public nuSQUIDS {
		private:
			squids::SU_vector NSI;
			std::vector<squids::SU_vector> NSI_evol;
			std::unique_ptr<double[]> hiBuffer;
			double HI_prefactor;

			// NSI parameters
			double epsilon_mutau;
			double epsilon_mumu;
			double epsilon_emu;
			double epsilon_ee;
			double epsilon_tautau;
			double epsilon_etau;

			void AddToPreDerive(double x);
			void AddToReadHDF5(hid_t hdf5_loc_id);
			void AddToWriteHDF5(hid_t hdf5_loc_id) const;

			squids::SU_vector HI(unsigned int ei,unsigned int index_rho) const;
			void Set_epsilons();
		public:
			//Constructor
			nuSQUIDSNSI(){;}
			//nuSQUIDSNSI(nuSQUIDSNSI&&);
			nuSQUIDSNSI(marray<double,1> Erange,unsigned int numneu, NeutrinoType NT=both, bool iinteraction=false, std::shared_ptr<NeutrinoCrossSections> ncs = nullptr);
			//nuSQUIDSNSI(unsigned int numneu, NeutrinoType NT = neutrino);
			//nuSQUIDSNSI(std::string hdf5_filename, std::string grp = "/",
			//		std::shared_ptr<InteractionStructure> int_struct = nullptr);
			//nuSQUIDSNSI& operator=(nuSQUIDSNSI&&);

			void print_matrix(const gsl_matrix_complex *m);
            // NSI parameter setters
			void Set_epsilon_mutau(double eps_mutau);
			void Set_epsilon_mumu(double eps_mumu);
			void Set_epsilon_tautau(double eps_tautau);
			void Set_epsilon_ee(double eps_ee);
			void Set_epsilon_emu(double eps_emu);
			void Set_epsilon_etau(double eps_etau);
			void Set_all_epsilons(double eps_mutau,double eps_mumu=0.0, double eps_tautau=0.0, double eps_emu=0.0,double eps_ee=0.0,double eps_etau=0.0);
	};
}
#endif
