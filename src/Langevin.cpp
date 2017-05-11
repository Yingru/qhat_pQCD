#include "Langevin.h"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include "qhat.h"

//giving vec(xi) in (px, py, pz) coordinate, Rotate vec(xi) to (px, py, py) --> (0, 0, p) coordinate
std::vector<double> rotate_toZ(std::vector<double> const& xi, double px, double py, double pz)
{
        double p_perp = std::sqrt(px*px + py*py);
        double p_length = std::sqrt(p_perp*p_perp + pz*pz);
        double cos_theta = pz / p_length;
        double sin_theta = p_perp / p_length;
        double cos_phi = px / p_perp;
        double sin_phi = py / p_perp;

        std::vector<double> Rxi;
        Rxi.resize(3);

        Rxi[0] = cos_theta*cos_phi*xi[0] + cos_theta*sin_phi*xi[1] - sin_theta*xi[2];
        Rxi[1] = - sin_phi*xi[0] + cos_phi * xi[1];
        Rxi[2] = sin_theta*cos_phi*xi[0] + sin_theta*sin_phi*xi[1] + cos_theta*xi[2];
        return Rxi;
}



//giving vec(xi) in (0, 0, p) coordinate, Rotate vec(xi) back to (px, py, pz) coordinate: (0, 0, p) --> (px, py, pz)
std::vector<double> rotate_fromZ(std::vector<double> const& xi, double px, double py, double pz)
{
        double p_perp = std::sqrt(px*px + py*py);

        if (p_perp < 1e-9) // in z direction already, no need for rotate
        {
                std::vector<double> Rxi(3);
                Rxi[0] = xi[0];
                Rxi[1] = xi[1];
                Rxi[2] = xi[2];
                return Rxi;
        }
        else // rotate
        {
                double p_length = std::sqrt(p_perp*p_perp + pz*pz);
                double cos_theta = pz / p_length;
                double sin_theta = p_perp / p_length;
                double cos_phi = px / p_perp;
                double sin_phi = py / p_perp;

                std::vector<double> Rxi(3);
                Rxi[0] = cos_theta*cos_phi*xi[0] - sin_phi*xi[1] + sin_theta*cos_phi*xi[2];
                Rxi[1] = cos_theta*sin_phi*xi[0] + cos_phi*xi[1] + sin_theta*sin_phi*xi[2];
                Rxi[2] = - sin_theta*xi[0] + cos_theta*xi[2];

                //std::cout << "Rxi: " << Rxi[0] << " " << Rxi[1] << " " << Rxi[2] << std::endl;
                return Rxi;
        }
}


int update_by_Langevin_test(particle& HQ, Qhat_2to2* qhatQq2Qq, Qhat_2to2* qhatQg2Qg, double temp, double deltat_lrf, bool EinR=true)
{
        double GeV_to_Invfm = 5.068;
        double drag_Qq, drag_Qg, drag, kperp_Qq, kperp_Qg, kperp, kpara_Qq, kpara_Qg, kpara;

        double p_perp = std::sqrt(HQ.p[1]*HQ.p[1] + HQ.p[2]*HQ.p[2]);
        double p_length = std::sqrt( p_perp*p_perp + HQ.p[3]*HQ.p[3]);
        double M2 = HQ.p[0]*HQ.p[0] - p_length*p_length;


        // step 1: pre_point shceme 
        double *args = new double[4];
        args[0] = HQ.p[0];
        args[1] = temp;
        args[2] = 0.0;
        args[3] = 1; drag_Qq = qhatQq2Qq->interpQ(args); drag_Qg = qhatQg2Qg->interpQ(args);
        args[3] = 2; kperp_Qq = qhatQq2Qq->interpQ(args); kperp_Qg = qhatQg2Qg->interpQ(args);
        args[3] = 3; kpara_Qq = qhatQq2Qq->interpQ(args); kpara_Qg = qhatQg2Qg->interpQ(args);
        delete [] args;

        drag = (drag_Qq + drag_Qg)/p_length * GeV_to_Invfm;
        kperp = (kperp_Qq + kperp_Qg) * GeV_to_Invfm;
        kpara = ((kpara_Qq - drag_Qq*drag_Qq) + (kpara_Qg - drag_Qg*drag_Qg) ) * GeV_to_Invfm;
        if (EinR)    drag = kperp/(2*temp*HQ.p[0]);
        //std::cout << "Langevin: " << HQ.p[0] << " " << drag << " " << kperp << " " << kpara << std::endl;

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);
        std::normal_distribution<double> distribution(0.0, 1.0);  // sample normal distribution with mean=0., width=1.0;

        std::vector<double> rho_xi(3, 0.0);
        for (size_t i=0; i<3; ++i)
                rho_xi[i] = distribution(generator);   // save the value to perform post_point scheme

        std::vector<double> xi(3, 0.0);
        xi[0] = std::sqrt(kperp/deltat_lrf) * rho_xi[0];
        xi[1] = std::sqrt(kperp/deltat_lrf) * rho_xi[1];
        xi[2] = std::sqrt(kpara/deltat_lrf) * rho_xi[2];


        std::vector<double> new_p(3, 0.0);  // remember, we always do the diffusion in (0,0,p_length) frame, and then rotate back
        new_p[0] = xi[0]*deltat_lrf; // x direction in (0,0,p_length) frame
        new_p[1] = xi[1]*deltat_lrf; // y direction in (0,0,p_length) frame
        new_p[2] = p_length + (-drag * p_length + xi[2])*deltat_lrf;


/*
        // step 2: post_point shceme
        double new_energy = std::sqrt(M2 + new_p[0]*new_p[0] + new_p[1]*new_p[1] + new_p[2]*new_p[2]);
        double new_drag_Qq, new_drag_Qg, new_drag, new_kperp_Qq, new_kperp_Qg, new_kperp, new_kpara_Qq, new_kpara_Qg, new_kpara;
        double *args_ = new double [4];

        args_[0] = new_energy;
        args_[1] = temp;
        args_[2] = 0.0;
        args_[3] = 1; new_drag_Qq = qhatQq2Qq->interpQ(args_); new_drag_Qg = qhatQg2Qg->interpQ(args_);
        args_[3] = 2; new_kperp_Qq = qhatQq2Qq->interpQ(args_); new_kperp_Qg = qhatQg2Qg->interpQ(args_);
        args_[3] = 3; new_kpara_Qq = qhatQq2Qq->interpQ(args_); new_kpara_Qg = qhatQg2Qg->interpQ(args_);
        delete [] args_;

        new_drag = (new_drag_Qq + new_drag_Qg) * GeV_to_Invfm;
        new_kperp = (new_kperp_Qq + new_kperp_Qg) * GeV_to_Invfm;
        new_kpara = ((new_kpara_Qq - new_drag_Qq*new_drag_Qq) + (new_kpara_Qg - new_drag_Qg*new_drag_Qg) ) * GeV_to_Invfm;

        xi[0] = std::sqrt(new_kperp/deltat_lrf) * rho_xi[0];
        xi[1] = std::sqrt(new_kperp/deltat_lrf) * rho_xi[1];
        xi[2] = std::sqrt(new_kpara/deltat_lrf) * rho_xi[2];
*/



        // now let us diffusion
        std::vector<double> p_inZ{0., 0., p_length};
        p_inZ[0] += xi[0] * deltat_lrf;
        p_inZ[1] += xi[1] * deltat_lrf;
        p_inZ[2] += (-drag*p_length + xi[2]) * deltat_lrf;

        std::vector<double> p_lrf(3);
        p_lrf = rotate_fromZ(p_inZ, HQ.p[1], HQ.p[2], HQ.p[3]);

        for (size_t i=0; i<3; ++i)
        {
                HQ.x[i] += HQ.p[i+1]/HQ.p[2] * deltat_lrf;
                HQ.p[i+1] = p_lrf[i];
        }

        HQ.p[0] = std::sqrt(M2 + HQ.p[1]*HQ.p[1] + HQ.p[2]*HQ.p[2] + HQ.p[3]*HQ.p[3]);

        return 0;
}




int update_by_Langevin(particle& HQ, double temp, double drag, double kpara, double kperp, double deltat_lrf, bool EinR=true)
// if EinR (Einstein relation is true... 
{
        //std::cout << "Langevin_old: " << HQ.p[0] << " " << drag << " " << kperp << " " << kpara << " " << std::endl;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator(seed);

        double p_perp = std::sqrt( HQ.p[1]*HQ.p[1] + HQ.p[2]*HQ.p[2]);
        double p_length = std::sqrt( p_perp*p_perp + HQ.p[3]*HQ.p[3]);
        double M2 = HQ.p[0]*HQ.p[0] - p_length*p_length;


        // here n (0, 0, p_length) coordinate
        std::vector<double> xi_width;
        xi_width.resize(3);
        xi_width[0] = std::sqrt(kperp/deltat_lrf);
        xi_width[1] = std::sqrt(kperp/deltat_lrf);
        xi_width[2] = std::sqrt(kpara/deltat_lrf);

        std::vector<double> xi_mean (3, 0.0);
        std::vector<double> xi(3, 0.0);

        for (int i=0; i<3; ++i)
        {
                std::normal_distribution<double> distribution(xi_mean[i], xi_width[i]);
                xi[i] = distribution(generator);
        }

        if (EinR)    drag = kperp / (2 * temp * HQ.p[0])  * p_length;  // here is \Gamma * pz
        
        std::vector<double> p_inZ {0., 0., p_length};
        p_inZ[0] += xi[0] * deltat_lrf;
        p_inZ[1] += xi[1] * deltat_lrf;
        p_inZ[2] += (-drag + xi[2]) * deltat_lrf;

        std::vector<double> p_lrf(3);
 
        //std::cout <<"p_inZ: "  <<  p_inZ[0] << " " << p_inZ[1] << " " << p_inZ[2] << std::endl;
        p_lrf = rotate_fromZ(p_inZ, HQ.p[1], HQ.p[2], HQ.p[3]);

        //std::cout << "p_lrf: " << p_lrf[0] << " " << p_lrf[1] << " " << p_lrf[2] << std::endl;
        for (int i=0; i<3; ++i)
        {
                HQ.p[i+1] = p_lrf[i];
                HQ.x[i] += HQ.p[i+1]/HQ.p[0] * deltat_lrf;
        }
        
        HQ.p[0] = std::sqrt(HQ.p[1]*HQ.p[1] + HQ.p[2]*HQ.p[2] + HQ.p[3]*HQ.p[3] + M2);
       
        //std::cout << HQ.p[0] << " " << HQ.p[1] << " " << HQ.p[2] << " " << HQ.p[3] << " " << M2 << std::endl;

        return 0;

}


/* although I am writing a large function here, actually, I can certainly make them piece by piece smaller*/


