#include <nuSQuIDS/nuSQuIDSNSI.h>

namespace nusquids{


    void nuSQUIDSNSI::AddToPreDerive(double x){
        for(int ei = 0; ei < ne; ei++){
            // assuming same hamiltonian for neutrinos/antineutrinos
            //SU_vector h0 = H0(E_range[ei],0);
            //NSI_evol[ei] = NSI.Evolve(h0,(x-Get_t_initial()));
            NSI_evol[ei] = NSI.Evolve(H0_array[ei],(x-Get_t_initial()));
        }
    }

    void nuSQUIDSNSI::AddToReadHDF5(hid_t hdf5_loc_id){
        // here we read the new parameters now saved in the HDF5 file
        hid_t nsi = H5Gopen(hdf5_loc_id, "nsi", H5P_DEFAULT);
        H5LTget_attribute_double(hdf5_loc_id,"nsi","mu_tau" ,&epsilon_mutau);
        H5LTget_attribute_double(hdf5_loc_id,"nsi","e_tau" ,&epsilon_etau);
        H5LTget_attribute_double(hdf5_loc_id,"nsi","e_mu" ,&epsilon_emu);
        H5LTget_attribute_double(hdf5_loc_id,"nsi","e_e" ,&epsilon_ee);
        H5LTget_attribute_double(hdf5_loc_id,"nsi","mu_mu" ,&epsilon_mumu);
        H5LTget_attribute_double(hdf5_loc_id,"nsi","tau_tau" ,&epsilon_tautau);
        H5Gclose(nsi);
    }

    void nuSQUIDSNSI::AddToWriteHDF5(hid_t hdf5_loc_id) const {
        // here we write the new parameters to be saved in the HDF5 file
        H5Gcreate(hdf5_loc_id, "nsi", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5LTset_attribute_double(hdf5_loc_id, "nsi","mu_tau",&epsilon_mutau, 1);
        H5LTset_attribute_double(hdf5_loc_id, "nsi","e_tau",&epsilon_etau, 1);
        H5LTset_attribute_double(hdf5_loc_id, "nsi","e_mu",&epsilon_emu, 1);
        H5LTset_attribute_double(hdf5_loc_id, "nsi","e_e",&epsilon_ee, 1);
        H5LTset_attribute_double(hdf5_loc_id, "nsi","mu_mu",&epsilon_mumu, 1);
        H5LTset_attribute_double(hdf5_loc_id, "nsi","tau_tau",&epsilon_tautau, 1);
    }

    squids::SU_vector nuSQUIDSNSI::HI(unsigned int ei, unsigned int index_rho) const{
        double CC = HI_prefactor*current_density*current_ye;

        // // construct potential in flavor basis
        squids::SU_vector potential(nsun,hiBuffer.get());

        potential = (3.0*CC)*NSI_evol[ei];

        if ((index_rho == 0 and NT==both) or NT==neutrino){
            // neutrino potential
            return nuSQUIDS::HI(ei,index_rho) + potential;
        } else if ((index_rho == 1 and NT==both) or NT==antineutrino){
            // antineutrino potential
            return nuSQUIDS::HI(ei,index_rho) + (-1.0)*std::move(potential);
        } else{
            throw std::runtime_error("nuSQUIDS::HI : unknown particle or antiparticle");
        }
    }

    //Constructor
    nuSQUIDSNSI::nuSQUIDSNSI(marray<double,1> Erange, unsigned int numneu, NeutrinoType NT, bool iinteraction,std::shared_ptr<NeutrinoCrossSections> ncs):
        nuSQUIDS(Erange, numneu, NT, iinteraction, ncs),
        hiBuffer(new double[nsun*nsun]), epsilon_mutau(0), epsilon_mumu(0), epsilon_tautau(0), epsilon_emu(0), epsilon_ee(0), epsilon_etau(0)
    {
        assert(numneu == 3);
        // defining a complex matrix M which will contain our flavor
        // violating flavor structure.
        // Added all epsilon flavor combinations
        Set_MixingParametersToDefault();
        Set_epsilons();
    }

    void nuSQUIDSNSI::print_matrix(const gsl_matrix_complex *m)
    {
        std::cout << std::endl;
        for (size_t i = 0; i < m->size1; i++) {
            for (size_t j = 0; j < m->size2; j++) {
                std::cout << m->data[(i * m->size2 + j)*2] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void nuSQUIDSNSI::Set_epsilon_mutau(double eps_mutau){
        epsilon_mutau = eps_mutau;
        Set_epsilons();
    }
    void nuSQUIDSNSI::Set_epsilon_mumu(double eps_mumu){
        epsilon_mumu = eps_mumu;
        Set_epsilons();
    }
    void nuSQUIDSNSI::Set_epsilon_tautau(double eps_tautau){
        epsilon_tautau = eps_tautau;
        Set_epsilons();
    }
    void nuSQUIDSNSI::Set_epsilon_ee(double eps_ee){
        epsilon_ee = eps_ee;
        Set_epsilons();
    }
    void nuSQUIDSNSI::Set_epsilon_emu(double eps_emu){
        epsilon_emu = eps_emu;
        Set_epsilons();
    }
    void nuSQUIDSNSI::Set_epsilon_etau(double eps_etau){
        epsilon_etau = eps_etau;
        Set_epsilons();
    }

    void nuSQUIDSNSI::Set_all_epsilons(double eps_mutau, double eps_mumu, double eps_tautau, double eps_emu, double eps_ee, double eps_etau){
        epsilon_mutau = eps_mutau;
        epsilon_mumu = eps_mumu;
        epsilon_tautau = eps_tautau;
        epsilon_ee = eps_ee;
        epsilon_emu = eps_emu;
        epsilon_etau = eps_etau;
        Set_epsilons();
    }

    void nuSQUIDSNSI::Set_epsilons(){
        gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3);
        gsl_complex c {{ epsilon_mutau , 0.0 }};
        gsl_complex cet {{ epsilon_etau , 0.0 }};
        gsl_complex cmm {{ epsilon_mumu , 0.0 }};
        gsl_complex ctt {{ epsilon_tautau , 0.0 }};
        gsl_complex cee {{ epsilon_ee , 0.0 }};
        gsl_complex cem {{ epsilon_emu , 0.0 }};
        gsl_matrix_complex_set(M,2,1,c);
        gsl_matrix_complex_set(M,1,2,gsl_complex_conjugate(c));
        gsl_matrix_complex_set(M,2,0,cet);
        gsl_matrix_complex_set(M,0,2,gsl_complex_conjugate(cet));
        gsl_matrix_complex_set(M,1,0,cem);
        gsl_matrix_complex_set(M,0,1,gsl_complex_conjugate(cem));
        gsl_matrix_complex_set(M,0,0,cee);
        gsl_matrix_complex_set(M,1,1,cmm);
        gsl_matrix_complex_set(M,2,2,ctt);
        NSI = squids::SU_vector(M);
        NSI.RotateToB1(params);
        NSI_evol.resize(ne);
        for(int ei = 0; ei < ne; ei++){
            NSI_evol[ei] = squids::SU_vector(nsun);
        }
        gsl_matrix_complex_free(M);
        HI_prefactor = params.sqrt2*params.GF*params.Na*pow(params.cm,-3);
    }
}
