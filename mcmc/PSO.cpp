#include "PSO.h"

PSO::PSO(manager *man) : LikelihoodFit(man) {

    fConstriction = fitMan->raw()["General"]["PSO"]["Constriction"].as<double>();
    fInertia = fitMan->raw()["General"]["PSO"]["Inertia"].as<double>()*fConstriction;
    fOne = fitMan->raw()["General"]["PSO"]["One"].as<double>()*fConstriction;
    fTwo = fitMan->raw()["General"]["PSO"]["Two"].as<double>()*fConstriction;
    fParticles = fitMan->raw()["General"]["PSO"]["Particles"].as<int>();
    fIterations = fitMan->raw()["General"]["PSO"]["Iterations"].as<int>();
    fConvergence = fitMan->raw()["General"]["PSO"]["Convergence"].as<double>();

    fDim = 0;

    if(fTestLikelihood)
    {
      fDim = fitMan->raw()["General"]["PSO"]["TestLikelihoodDim"].as<int>();
    }
}

void PSO::runMCMC(){

    PrepareFit();

    if(fTestLikelihood){ 
        outTree->Branch("nParts", &fParticles, "nParts/I");
        for(int i = 0; i < fDim; ++i){
            par = new double[fParticles];
            paramlist.push_back(par);
            outTree->Branch(Form("Parameter_%d", i), paramlist[i], Form("Parameter_%d[nParts]/D",i));
        }
//        vel = new double[fParticles];
        outTree->Branch("vel", vel, "vel[nParts]/D");
    }

    init();
    run();
    WriteOutput();
    return;
}


void PSO::init(){

    fBestValue = 1234567890.0;

    //KS: For none PCA this will be eqaul to normal parameters
    //const int NparsPSOFull = NPars;
    //const int NparsPSO = NParsPCA;

    std::cout << "Preparing PSO" << std::endl;

    // Initialise bounds on parameters
    if(fTestLikelihood){
        for (int i = 0; i < fDim; i++){
            // Test function ranges
            ranges_min.push_back(-5);
            ranges_max.push_back(5);
            fixed.push_back(0);
        }
    }
    else{
        for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it){
            if(!(*it)->IsPCA())
            {
                fDim += (*it)->GetNumParams();
                for(int i = 0; i < (*it)->GetNumParams(); ++i)
                {
                    double curr = (*it)->getParInit(i);
                    double lim = 10.0*(*it)->getDiagonalError(i);
                    double low = (*it)->GetLowerBound(i);
                    double high = (*it)->GetUpperBound(i);
                    if(low > curr - lim) ranges_min.push_back(low);
                    else ranges_min.push_back(curr - lim);
                    if(high < curr + lim) ranges_min.push_back(high);
                    else ranges_min.push_back(curr + lim);
                    prior.push_back(curr);

                    if((*it)->isParameterFixed(i)){
                        fixed.push_back(1);
                    }
                    else{
                        fixed.push_back(0);
                    }
                }
            }
            else
            {
                fDim += (*it)->getNpars();
                for(int i = 0; i < (*it)->getNpars(); ++i)
                {
                    ranges_min.push_back(-100.0);
                    ranges_max.push_back(100.0);
                    prior.push_back((*it)->getParInit(i));
                    if((*it)->isParameterFixedPCA(i)){
                        fixed.push_back(1);
                    }
                    else{
                        fixed.push_back(0);
                    }
                }
            }
        }
    }

    std::cout << "Printing Minimums and Maximums of Variables to be minimised" << std::endl;
    for (int i =0; i<fDim; i++){
        std::cout << "Variable "<< i<<" :" << ranges_min[i] << ", "<< ranges_max[i] << std::endl;
    }

    // Initialise particle positions
    for (int i = 0; i < fParticles; ++i){
        std::vector<double> init_position;
        std::vector<double> init_velocity;

        //Initialising in +/- 5sigma of prior value from BANFF interface
        for (int j=0; j<fDim; ++j){
            if(fixed[j]){
                init_position.push_back(prior[j]);
                init_velocity.push_back(0.0);
            }
            else{
                double dist = fabs(ranges_max[j]-ranges_min[j]);
                //Initialise to random position uniform in space
                init_position.push_back(ranges_min[j] + random->Rndm()*dist);
                //Initialise velocity to random position uniform in space
                init_velocity.push_back((2.0*random->Rndm()-1.0));//*dist);
            }
        }

        particle* new_particle = new particle(init_position,init_velocity);
        new_particle->set_personal_best_position(init_position);
        double new_value = CalcChi(init_position);
        new_particle->set_personal_best_value(new_value);
        new_particle->set_value(new_value);
        system.push_back(new_particle);
        if(new_value < fBestValue){
            fBestValue = new_value;
            set_best_particle(new_particle);
        }
    }
}

std::vector<std::vector<double> > PSO::bisection(std::vector<double>position,double minimum, double range, double precision){
    std::vector<std::vector<double>> uncertainties_list;
    for (unsigned int i = 0; i< position.size(); ++i){
        std::cout << i << std::endl;
        //std::vector<double> uncertainties;
        std::vector<double> new_position = position; new_position[i] = position[i]-range;
        double val_1 = CalcChi(new_position)-minimum-1.0;
        while (val_1*-1.0> 0.0){
            new_position[i] -= range;
            val_1 = CalcChi(new_position)-minimum-1.0;
        }
        std::vector<double> bisect_position = position; bisect_position[i] = bisect_position[i] - (position[i]-new_position[i])/2;
        std::vector<std::vector<double>> position_list{new_position,bisect_position,position};
        double val_2 = CalcChi(bisect_position)-minimum-1.0;
        std::vector<double> value_list{val_1,val_2, -1.0};
        double res = 1.0;
        while (res > precision){
            if (value_list[0] * value_list[1] < 0){
                std::vector<double> new_bisect_position = position_list[0];new_bisect_position[i] =new_bisect_position[i]+ (position_list[1][i]-position_list[0][i])/2;
                double new_val = CalcChi(new_bisect_position)-minimum-1.0;
                position_list[2] = position_list[1];
                value_list[2] = value_list[1];
                value_list[1] = new_val;
                position_list[1] = new_bisect_position;
                res = abs(position[2]-position[0]);
            }
            else{
                std::vector<double> new_bisect_position = position_list[1];new_bisect_position[i] += (position_list[2][i]-position_list[1][i])/2;
                double new_val = CalcChi(new_bisect_position)-minimum-1.0;
                position_list[0] = position_list[1];
                value_list[0] = value_list[1];
                value_list[1] = new_val;
                position_list[1] = new_bisect_position;
                res = abs(position_list[2][i]-position_list[1][i]);
            }
        }
        //do the same thing for position uncertainty
        std::vector<double> new_position_p = position; new_position_p[i] = position[i]+range;
        double val_1_p = CalcChi(new_position_p)-minimum-1.0;
        while (val_1_p * -1.0 > 0.0){
            new_position_p[i] += range;
            val_1_p = CalcChi(new_position_p)-minimum-1.0;
        }
        std::vector<double> bisect_position_p = position; bisect_position_p[i] = bisect_position_p[i] += (new_position_p[i]-position[i])/2;
        std::vector<std::vector<double>> position_list_p{position,bisect_position_p,new_position_p};
        double val_2_p = CalcChi(bisect_position_p)-minimum-1.0;
        std::vector<double> value_list_p{-1.0,val_2_p, val_1_p};
        double res_p = 1.0;
        while (res_p > precision){
            if (value_list_p[0] * value_list_p[1] < 0){
                std::vector<double> new_bisect_position_p = position_list_p[0];new_bisect_position_p[i] += (position_list_p[1][i]-position_list_p[0][i])/2;
                double new_val_p = CalcChi(new_bisect_position_p)-minimum-1.0;
                position_list_p[2] = position_list_p[1];
                value_list_p[2] = value_list_p[1];
                value_list_p[1] = new_val_p;
                position_list_p[1] = new_bisect_position_p;
                res = abs(position[2]-position[0]);
                res_p = abs(position_list_p[1][i]-position_list_p[0][i]);
                //std::cout << "Pos midpoint is " << position_list_p[1][i] << std::endl;
            }
            else{
                std::vector<double> new_bisect_position_p = position_list_p[1];new_bisect_position_p[i] += (position_list_p[2][i]-position_list_p[1][i])/2;
                double new_val_p = CalcChi(new_bisect_position_p)-minimum-1.0;
                position_list_p[0] = position_list_p[1];
                value_list_p[0] = value_list_p[1];
                value_list_p[1] = new_val_p;
                position_list_p[1] = new_bisect_position_p;
                res_p = abs(position_list_p[2][i]-position_list_p[1][i]);
                //std::cout << "Pos midpoint is " << position_list_p[1][i] << std::endl;
            }
        }
        uncertainties_list.push_back({abs(position[i]-position_list[1][i]),abs(position[i]-position_list_p[1][i])});
        std::cout << "Uncertainty finished for d = "<< i << std::endl;
        std::cout << std::setprecision(10)<< "LLR values for ± positive and negative uncertainties are " << CalcChi(position_list[1]) << " and " << CalcChi(position_list_p[1]) << std::endl;
    }
    return uncertainties_list;
}

std::vector<std::vector<double>> PSO::calc_uncertainty(std::vector<double>position,double minimum) {
    std::vector<double> pos_uncertainty(position.size());
    std::vector<double> neg_uncertainty(position.size());
    int num = 200;
    std::vector<double> pos = position;
    for (unsigned int i = 0; i < position.size(); ++i) {

        double curr_ival = pos[i];

        double neg_stop = position[i] - 5e-2;
        double pos_stop = position[i] + 5e-2;
        double start = position[i];
        std::vector<double> x(num);
        std::vector<double> y(num);
        double StepPoint = (start-neg_stop) / (num - 1);
        double value = start;
        for (int j = 0; j < num; ++j) {
          pos[i] = value;
          double LLR = CalcChi(position) - minimum - 1.0;
          x[j] = value;
          y[j] = LLR;
          value -= StepPoint;
        }
        pos[i] = curr_ival;

        int closest_index = 0;
        double closest_value = abs(y[0]); // Initialize with the first element
        for (unsigned int ii = 1; ii < y.size(); ++ii) {
          double abs_y = abs(y[ii]);
          if (abs_y < closest_value) {
            closest_index = ii;
            closest_value = abs_y;
          }
        }
        neg_uncertainty[i] = x[closest_index];
        std::cout << "Neg" << std::endl;
        x.assign(num, 0);
        y.assign(num, 0);
        StepPoint = (pos_stop-start) / (num - 1);
        value = start;
        for (int j = 0; j < num; ++j) {
          pos[i] = value;
          double LLR = CalcChi(position) - minimum - 1.0;
          x[j] = value;
          y[j] = LLR;
          value += StepPoint;
        }
        pos[i] = curr_ival;
        closest_index = 0;
        closest_value = abs(y[0]); // Initialize with the first element
        for (unsigned int ii = 1; ii < y.size(); ++ii) {
          double abs_y = abs(y[ii]);
          if (abs_y < closest_value) {
            closest_index = ii;
            closest_value = abs_y;
          }
        } 
        pos_uncertainty[i] = x[closest_index];
    } 
    std::vector<std::vector<double>> res{neg_uncertainty,pos_uncertainty};
    return res;
}

void PSO::uncertainty_check(std::vector<double> previous_pos){
    std::vector<std::vector<double >> x_list;
    std::vector<std::vector<double >> y_list;
    std::vector<double> position = previous_pos;
    int num = 5000;
    for (unsigned int i = 0;i<previous_pos.size();++i){
        double curr_ival = position[i];
        double start = previous_pos[i] - 1e-1;
        double stop = previous_pos[i] + 1e-1;
        std::vector<double> x(num);
        std::vector<double> y(num);
        double StepPoint = (stop - start) / (num - 1);
        double value = start;
        // std::cout << "result for fDim " << 1 << std::endl;
        for (int j =0;j< num; ++j){
            position[i] = value;
            double LLR = CalcChi(position);
            x[j] = value;
            y[j] = LLR;
            value += StepPoint;
        }
        position[i] = curr_ival;
        std::cout << " " << std::endl;
        std::cout << "For fDim" << i+1  << " x list is " ;
        for (unsigned int k= 0;k<x.size(); ++k){
            std::cout << x[k] << " , " ;
        } std::cout << " " << std::endl;
        std::cout << "  " << std::endl;
        std::cout << "For fDim" << i+1 << " y list is " ;
        for (unsigned int k= 0;k<x.size(); ++k){
            std::cout << y[k] << " , " ;
        } std::cout <<  " " << std::endl;
        std::cout << " " <<std::endl;
    }
}

double PSO::swarmIterate(){

    std::vector<double> total_pos(fDim,0.0);

    for (int i = 0; i < fParticles; ++i) {

        std::vector<double> part1 = vector_multiply(system[i]->get_velocity(), fInertia);
        std::vector<double> part2 = vector_multiply(vector_subtract(system[i]->get_personal_best_position(), system[i]->get_position()), (fOne * random->Rndm()));
        std::vector<double> part3 = vector_multiply(vector_subtract(get_best_particle()->get_personal_best_position(), system[i]->get_position()),(fTwo * random->Rndm()));
        std::vector<double> new_velocity = three_vector_addition(part1, part2, part3);
        std::vector<double> new_pos = vector_add(system[i]->get_position(), new_velocity);
        transform(total_pos.begin(), total_pos.end(), new_pos.begin(), total_pos.begin(),[](double x, double y) {return x+y;});

        for (int j = 0; j < fDim; ++j) {
            // Check if out of bounds and reflect if so
            if(ranges_min[j] > new_pos[j]){
                new_pos[j] = ranges_min[j];
            }
            else if(new_pos[j] > ranges_max[j]) {
                new_pos[j] = ranges_max[j];
            }
            // If parameter fixed don't update it
            if(fixed[j]) new_pos[j] = system[i]->get_position()[j];
        }

        if(fTestLikelihood){
            double velo = 0.0;
            for (int j = 0; j < fDim; ++j) {
                paramlist[j][i] = new_pos[j];
                velo += new_velocity[j]*new_velocity[j];
            }
            vel[i] = sqrt(velo);
        }

        system[i]->set_velocity(new_velocity);
        system[i]->set_position(new_pos);
        double new_value = CalcChi(new_pos);
        if(new_value <= system[i]->get_personal_best_value()) {
            system[i]->set_personal_best_value(new_value);
            system[i]->set_personal_best_position(new_pos);
            if(new_value < fBestValue){
                fBestValue = new_value;
                set_best_particle(system[i]);
            }
        }
    }

    std::vector<double> best_pos = get_best_particle()->get_personal_best_position();
    std::vector<double> result(best_pos.size(), 0.0);
    transform(total_pos.begin(), total_pos.end(), total_pos.begin(), [=](double val){return val/fParticles;});
    transform(total_pos.begin(),total_pos.end(),best_pos.begin(),result.begin(),[](double x, double y) {return x-y;});

    double mean_dist_sq = 0;
    for (int i = 0; i<fDim; i++){
        mean_dist_sq += result[i]*result[i];
    }

    return mean_dist_sq;
}

void PSO::run() {

    double mean_dist_sq = 0;

    int iter = 0;
    for(int i = 0; i < fIterations; ++i, ++iter){
        mean_dist_sq = swarmIterate();
        //double meanVel = std::accumulate(vel, vel + fParticles, 0) / fParticles;

        // Weight inertia randomly but scaled by total distance of swarm from global minimum - proxy for total velocity
        // fWeight = ((random->Rndm()+1.0)*0.5)*(10.0/meanVel);

        logLCurr = fBestValue;

        outTree->Fill();
        // Auto save the output
        if (step % auto_save == 0) outTree->AutoSave();
        step++;
        accCount++;

        if (i%100 == 0){
            std::cout << "Mean Dist Sq = " << mean_dist_sq <<std::endl;
            std::cout << "Current LLR = " << fBestValue << std::endl;
            std::cout << "Position = " <<std::endl;
            for (int j = 0; j< fDim; ++j){
                std::cout << "    Dim " << j << " = " << std::setprecision(10) << get_best_particle()->get_personal_best_position()[j] << std::endl;
            }
            
        }
        if(fConvergence > 0.0){
            if(mean_dist_sq < fConvergence){
                break;
            }
        }
    }
    std::cout << "Finished after " << iter <<" runs out of "<< fIterations << std::endl;
    std::cout << "Mean Dist " << mean_dist_sq <<std::endl;
    std::cout << "Best LLR " << get_best_particle()->get_personal_best_value() << std::endl;

    uncertainties = bisection(get_best_particle()->get_personal_best_position(),get_best_particle()->get_personal_best_value(),0.5,0.005);
    std::cout << "Position for Global Minimum = "<<std::endl;
    for (int i = 0; i< fDim; ++i){
        std::cout << "    Dim " << i << " = " << std::setprecision(10) << get_best_particle()->get_personal_best_position()[i]  << " +" << uncertainties[i][1] << ", -" <<  uncertainties[i][0] << std::endl;
    }
}

void PSO::WriteOutput(){

    outputFile->cd();

    TVectorD* PSOParValue = new TVectorD(fDim);
    TVectorD* PSOParError = new TVectorD(fDim);

    for(int i = 0; i < fDim; ++i)
    {
        (*PSOParValue)(i) = 0;
        (*PSOParError)(i) = 0;
    }

    std::vector<double> minimum = get_best_particle()->get_personal_best_position();

    int ParCounter = 0;

    if(fTestLikelihood){
        for(int i = 0; i < fDim; ++i){
            (*PSOParValue)(i) = minimum[i];
            (*PSOParError)(i) = (uncertainties[i][0]+uncertainties[i][1])/2.0;
        }
    }
    else{
        for (std::vector<covarianceBase*>::iterator it = systematics.begin(); it != systematics.end(); ++it)
        {
            if(!(*it)->IsPCA())
            {
                for(int i = 0; i < (*it)->GetNumParams(); ++i, ++ParCounter)
                {
                    double ParVal = minimum[ParCounter];
                    //KS: Basically apply mirroring for parameters out of bounds
                    (*PSOParValue)(ParCounter) = ParVal;
                    (*PSOParError)(ParCounter) = (uncertainties[ParCounter][0]+uncertainties[ParCounter][1])/2.0;
                    //KS: For fixed params HESS will not calcuate error so we need to pass prior error
                    if((*it)->isParameterFixed(i))
                    {
                        (*PSOParError)(ParCounter) = (*it)->getDiagonalError(i);
                    }
                }
            }
            else
            {
                //KS: We need to convert parameters from PCA to normal base
                TVectorD ParVals((*it)->GetNumParams());
                TVectorD ParVals_PCA((*it)->getNpars());

                TVectorD ErrorVals((*it)->GetNumParams());
                TVectorD ErrorVals_PCA((*it)->getNpars());

                //First save them
                //KS: This code is super convoluted as MaCh3 can store separate matrices while PSO has one matrix. In future this will be simplified, keep it like this for now.
                const int StartVal = ParCounter;
                for(int i = 0; i < (*it)->getNpars(); ++i, ++ParCounter)
                {
                    ParVals_PCA(i) = minimum[ParCounter];
                    ErrorVals_PCA(i) = (uncertainties[ParCounter][0]+uncertainties[ParCounter][1])/2.0;
                }
                ParVals = ((*it)->getTransferMatrix())*ParVals_PCA;
                ErrorVals = ((*it)->getTransferMatrix())*ErrorVals_PCA;

                ParCounter = StartVal;
                //KS: Now after going from PCA to normal let';s save it
                for(int i = 0; i < (*it)->GetNumParams(); ++i, ++ParCounter)
                {
                    (*PSOParValue)(ParCounter) = ParVals(i);
                    (*PSOParError)(ParCounter) = std::fabs(ErrorVals(i));
                    //int ParCounterMatrix = StartVal;
                    //If fixed take prior
                    if((*it)->isParameterFixedPCA(i))
                    {
                        (*PSOParError)(ParCounter) = (*it)->getDiagonalError(i);
                    }
                }
            }
        }
    }

    PSOParValue->Write("PSOParValue");
    PSOParError->Write("PSOParError");
    delete PSOParValue;
    delete PSOParError;
    // Save all the output
    SaveOutput();
}



// *******************
double PSO::CalcChi2(const double* x) {
// *******************

  if(fTestLikelihood)
  {
    return rastriginFunc(x);
  }
  else
  {
    return LikelihoodFit::CalcChi2(x);
  }
}


// *************************
double PSO::rastriginFunc(const double* x) {
// *************************

  stepClock->Start();

  //Search range: [-5.12, 5.12]
  const double A = 10.0;
  double sum = 0.0;
  for (int i = 0; i < fDim; ++i) {
      sum += x[i] * x[i] - A * cos(2.0 * 3.14 * x[i]);
  }
  double llh = A * fDim + sum;

  accProb = 1;

  stepClock->Stop();
  stepTime = stepClock->RealTime();

  return llh;
}
