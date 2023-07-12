//////////////////////////////////////////////////////////////////////////
////////////////        alpine_rider.cxx             /////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example            /////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:         Alp rider problem                 ////////////////
//////// Last modified: 09 February 2009                  ////////////////
//////// Reference:     Betts (2001)                      ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which////////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
//#define pi 4.0*atan(1.0)
/////////////////////////////////////////////aer) cost function //////////
//////////////////////////////////////////////////////////////////////////
adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return 0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{

    adouble dx1= states[4];
    adouble dx2= states[5];
    adouble dx3= states[6];
    adouble dx4= states[7];
    adouble u1 = controls[0];
    adouble u2 = controls[1];

    adouble L;

    // Eigen::Matrix<adouble,6,4> P;
    //     P <<-1.35,1,-0.5,0,
    //          1,0,0,0,
    //          0,1,0,0,
    //          0,0,1,0,
    //          0,1,0.5,1.35,
    //          0,0,0,1;
    // Eigen::Matrix<adouble,6,2> Fp;
    //     Fp <<0,0,
    //         1,0,
    //         0,0,
    //         -1,-1,
    //         0,0,
    //         0,1;
    L = -u1*(dx3-dx1)-u2*(dx4-dx3);

    return  L;
}

// adouble func2(adouble L){
//     return L;
// }

//求定积分//
// adouble integral(adouble(*f)(adouble), adouble min, adouble max){
// 	adouble result = 0.0;
// 	const int N = 1000;
// 	adouble delta = (max - min) / N;
//     Eigen::MatrixXf J(1, N);
// 	for (adouble i = min+delta; i < max; i+=delta)
// 	{
//         J=result;
// 		result += f(i)*delta;
// 	}
// 	return result;
// }

///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    Eigen::Map<Eigen::Matrix<adouble,8,1>> x(states,8);
    Eigen::Map<Eigen::Matrix<adouble,2,1>> u(controls,2);
    Eigen::Map<Eigen::Matrix<adouble,8,1>> dx(derivatives,8);
    adouble x1= states[0];
    adouble x2= states[1];
    adouble x3= states[2];
    adouble x4= states[3];
    adouble dx1=states[4];
    adouble dx2=states[5];
    adouble dx3=states[6];
    adouble dx4=states[7];
    adouble u1 = controls[0];
    adouble u2 = controls[1];
    adouble t  = time;

    //Define the model matrices A and B
    Eigen::Matrix<adouble,8,8> A;
        A <<0,0,0,0,1.0,1.0,1.0,1.0,
            0,0,0,0,1.0,1.0,1.0,1.0,
            0,0,0,0,1.0,1.0,1.0,1.0,
            0,0,0,0,1.0,1.0,1.0,1.0,
            -19.28,10.78,-4.54,0.15,-0.32,0.19,-0.12,-0.008,
            -0.33,-4.57,-0.001,0.33,0.017,-0.08,0,-0.02,
            0.23,0,-7.0,0.23,-0.012,0,-0.0096,-0.012,
            0.15,-10.78,-4.54,-19.28,-0.008,-0.188,-0.12,-0.3;
    Eigen::Matrix<adouble,8,4> B;
        B <<0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            0,0,0,0,
            0.0015,0.0007,-0.0005,-0.0003,
            0.0007,0.001,0,-0.0007,
            -0.0005,0,0.0014,-0.0005,
            -0.0003,-0.0007,-0.0005,0.0015;
    Eigen::Matrix<adouble,4,1> Fe;
        Fe << -7163.5*sin(4.0*atan(1.0)*0.5*t),
                13004.5*sin(4.0*atan(1.0)*0.5*t),
                1333.5*sin(4.0*atan(1.0)*0.5*t),
                9152*sin(4.0*atan(1.0)*0.5*t);
    Eigen::Matrix<adouble,4,1> Fpto;
        Fpto<<-u1,
                0,
                u1-u2,
                u2;
    dx=A*x+B*(Fe+Fpto);


}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x1i = initial_states[ 0 ];
   adouble x2i = initial_states[ 1 ];
   adouble x3i = initial_states[ 2 ];
   adouble x4i = initial_states[ 3 ];
   adouble x5i = initial_states[ 4 ];
   adouble x6i = initial_states[ 5 ];
   adouble x7i = initial_states[ 6 ];
   adouble x8i = initial_states[ 7 ];



   e[ 0 ] = x1i;
   e[ 1 ] = x2i;
   e[ 2 ] = x3i;
   e[ 3 ] = x4i;
   e[ 4 ] = x5i;
   e[ 5 ] = x6i;
   e[ 6 ] = x7i;
   e[ 7 ] = x8i;

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

    // adouble result=0;
    // result=integral(func2,0,50.0);
    // cout << result << endl;

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Alp rider problem";
    problem.outfilename                 = "alpine.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 8;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 8;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes                     << 240;



    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.phases(1).bounds.lower.states <<    -4, -4, -4, -4, -4, -4, -4, -4;

    problem.phases(1).bounds.upper.states <<     4, 4, 4, 4, 4, 4, 4, 4;

    problem.phases(1).bounds.lower.controls <<  -600.0, -300 ;

    problem.phases(1).bounds.upper.controls <<   600.0,  300 ;

    problem.phases(1).bounds.lower.events <<  0, 0, 0, 0, 0, 0, 0, 0;

    problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;

    // problem.phases(1).bounds.upper.path         << 100.0;

    // problem.phases(1).bounds.lower.path         << 0.0;

    problem.phases(1).bounds.lower.StartTime    = 0.0;

    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 20.0;

    problem.phases(1).bounds.upper.EndTime      = 20.0;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	   = &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		      = &events;
    problem.linkages		      = &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			               = problem.phases(1).nodes(0);

    MatrixXd x_guess    		            =  zeros(8,nnodes);

    x_guess.row(0) 			               = linspace(1,-1,nnodes);
    x_guess.row(1) 			               = linspace(1,-1,nnodes);
    x_guess.row(2) 			               = linspace(1,-1,nnodes);
    x_guess.row(3) 			               = linspace(1,-1,nnodes);
    x_guess.row(4)                         = linspace(1,-1,nnodes);
    x_guess.row(5)                         = linspace(1,-1,nnodes);
    x_guess.row(6)                         = linspace(1,-1,nnodes);
    x_guess.row(7)                         = linspace(1,-1,nnodes);

    problem.phases(1).guess.controls      = zeros(2,nnodes);
    problem.phases(1).guess.states        = x_guess;
    problem.phases(1).guess.time          = linspace(0.0,20.0,nnodes+1);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    // algorithm.nlp_iter_max                = 1000;
    // algorithm.nlp_tolerance               = 1.e-6;
    // algorithm.nlp_method                  = "IPOPT";
    // algorithm.scaling                     = "automatic";
    // algorithm.derivatives                 = "automatic";
    // algorithm.jac_sparsity_ratio          = 0.20;
    // algorithm.collocation_method          = "Legendre";
    // algorithm.diff_matrix                 = "central-differences";
    // algorithm.mesh_refinement             = "automatic";
    // algorithm.mr_max_increment_factor     = 0.3;
    // algorithm.mr_max_iterations           = 3;
    // algorithm.defect_scaling              = "jacobian-based";

    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.nlp_iter_max                = 5000;
    algorithm.nlp_tolerance               = 1.e-6;




////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   //////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    MatrixXd x_ = solution.get_states_in_phase(1);
    MatrixXd u_ = solution.get_controls_in_phase(1);
    MatrixXd t_ = solution.get_time_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

   Save(x_,"x.dat");
   Save(u_,"u.dat");
   Save(t_,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    plot(t_,x_.row(0),problem.name+": state", "time (s)", "state","x1");

    plot(t_,x_.row(1),problem.name+": state", "time (s)", "state","x2");

    plot(t_,x_.row(2),problem.name+": state", "time (s)", "state","x3");

    plot(t_,x_.row(3),problem.name+": state", "time (s)", "state","x4");

    plot(t_,x_.row(4),problem.name+": state", "time (s)", "state","x5");

    plot(t_,x_.row(5),problem.name+": state", "time (s)", "state","x6");

    plot(t_,x_.row(6),problem.name+": state", "time (s)", "state","x7");

    plot(t_,x_.row(7),problem.name+": state", "time (s)", "state","x8");

    plot(t_,x_.row(7)-x_.row(6),problem.name+": state", "time (s)", "state","theta2");

    plot(t_,x_.row(6)-x_.row(4),problem.name+": state", "time (s)", "state","theta1");


    plot(t_,u_.row(0),problem.name+": control","time (s)", "control", "u1");

    plot(t_,u_.row(1),problem.name+": control","time (s)", "control", "u2");

    plot(t_,x_.row(0),problem.name+": state x1", "time (s)", "state","x1",
                                        "pdf", "alpine_state1.pdf");

    plot(t_,x_.row(1),problem.name+": state x2", "time (s)", "state","x2",
                                        "pdf", "alpine_state2.pdf");

    plot(t_,x_.row(2),problem.name+": state x3", "time (s)", "state","x3",
                                        "pdf", "alpine_state3.pdf");

    plot(t_,x_.row(3),problem.name+": state x4", "time (s)", "state","x4",
                                        "pdf", "alpine_state4.pdf");
    
    plot(t_,x_.row(4),problem.name+": state x5", "time (s)", "state","x5",
                                        "pdf", "alpine_state5.pdf");

    plot(t_,x_.row(5),problem.name+": state x6", "time (s)", "state","x6",
                                        "pdf", "alpine_state6.pdf");

    plot(t_,x_.row(6),problem.name+": state x7", "time (s)", "state","x7",
                                        "pdf", "alpine_state7.pdf");

    plot(t_,x_.row(7),problem.name+": state x8", "time (s)", "state","x8",
                                        "pdf", "alpine_state8.pdf");                                    
                                        
    plot(t_,u_.row(0),problem.name+": control u1","time (s)", "control", "u1",
                                        "pdf", "alpine_control1.pdf");

    plot(t_,u_.row(1),problem.name+": control u1","time (s)", "control", "u2",
                                        "pdf", "alpine_control2.pdf");

    plot(t_,x_.row(7)-x_.row(6),problem.name+": state theta2", "time (s)", "theta2","x8",
                                        "pdf", "alpine_state9.pdf");

    plot(t_,x_.row(6)-x_.row(4),problem.name+": state theta1", "time (s)", "theta1","x8",
                                        "pdf", "alpine_state10.pdf");

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
