#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "rebound.h"


void heartbeat(struct reb_simulation* r);
void basic_save(struct reb_simulation* r, char* filename_a, char* filename_e, char* filename_inc, char* filename_omega, char* filename_Omega, char* filename_M, char* filename1_a, char* filename1_e, char* filename1_inc, char* filename1_omega, char* filename1_Omega, char* filename1_M);

// sim settings
int orbits = 1e5; //total number of orbits
double orbit_time = 2*M_PI;
double tmax;
double disk_mass = 1e-2; // Total disk mass
double bowl_mass = 0; // Total bowl mass
double BH_mass = 1.; // one solar mass
int num_bin_saves = 10; // number of binary simulation saves per sim
int num_saves_per_orbit = 10; // the number of times per orbit that the orbital elements get saved
int num_disk = 10; // total number of disk particles
int num_bowl = 10; // total number of bowl particles
int bin_save_count;
double orbit_counter = 0.;
double save_counter = 0.;
char tfilename[100];

char ofilename_a[100];
char ofilename_e[100];
char ofilename_inc[100];
char ofilename_omega[100];
char ofilename_Omega[100];
char ofilename_M[100];

char ofilename1_a[100];
char ofilename1_e[100];
char ofilename1_inc[100];
char ofilename1_omega[100];
char ofilename1_Omega[100];
char ofilename1_M[100];

char efilename[100];
char binprefix[100];
int bowl_type = 0; // 0 no Bowl, 1 bowl of black holes, 2 bowl of stars simular to milky way
clock_t last_hb;

int main(int argc, char* argv[]){
    // Create simulation
	struct reb_simulation* r = reb_create_simulation();
    /////////////////////////////////////////////////////////////////////////

	// Setup constants
    r->G			    = 1;		    // Gravitational constant, sets the units
    r->integrator		= REB_INTEGRATOR_IAS15;
    r->heartbeat		= heartbeat;
    
	//Disk setup
    double part_mD = disk_mass / (double)num_disk;        //All have equal mass
    double part_aD[num_disk];
    double part_eD[num_disk];
    double part_iD[num_disk];
    double part_oD[num_disk];
    double part_OD[num_disk];
    double part_MD[num_disk];
    for (int p = 0; p < num_disk; p++) {
        part_aD[p] = reb_random_uniform(.9, 1.9); //inverse square ditribution in a
        part_eD[p] = 0.7;                             //Identical eccentricity
        part_iD[p] = 1e-4;
        part_oD[p] = 0; //evenly distributed in argument of periapsis
        part_OD[p] = 0; //evenly distributed in longitude of periapsis
        part_MD[p] = reb_random_uniform(0., 2.*M_PI); //evenly distributed in mean anomaly
    }

	// Bowl setup
	double part_mB = bowl_mass / (double)num_bowl;        //All have equal mass
	double part_aB[num_bowl];
	double part_eB[num_bowl];
	double part_iB[num_bowl];
	double part_oB[num_bowl];
	double part_OB[num_bowl];
	double part_MB[num_bowl];
	if (bowl_type != 0){
		if (bowl_type == 1){
			for (int p = 0; p < num_bowl; p++) {
				part_aB[p] = reb_random_uniform(.2, .3); //inverse square ditribution in a
				part_eB[p] = 0.7;                      	//eccentricity
				part_iB[p] = reb_random_uniform(0., 2.*M_PI); //inclination
				part_oB[p] = reb_random_uniform(0., 2.*M_PI); //evenly distributed in argument of periapsis
				part_OB[p] = reb_random_uniform(0., 2.*M_PI); //evenly distributed in longitude of periapsis
				part_MB[p] = reb_random_uniform(0., 2.*M_PI); //evenly distributed in mean anomaly
			}
		}
		if (bowl_type == 2){
			for (int p = 0; p < num_bowl; p++) {
				part_aB[p] = reb_random_uniform(.01, 1.); //inverse square ditribution in a
				part_eB[p] = reb_random_uniform(0.3, 0.9);  //eccentricity
				part_iB[p] = reb_random_uniform(0., 2.*M_PI); //inclination
				part_oB[p] = reb_random_uniform(0., 2.*M_PI); //evenly distributed in argument of periapsis
				part_OB[p] = reb_random_uniform(0., 2.*M_PI); //evenly distributed in longitude of periapsis
				part_MB[p] = reb_random_uniform(0., 2.*M_PI); //evenly distributed in mean anomaly
			}
		}

	}
	/////////////////////////////////////////////////////////////////////////
    // Create the star
	struct reb_particle BH    = {0}; 
	BH.m                      = BH_mass;
    BH.hash                   = reb_hash("BH");
	reb_add(r, BH);
    /////////////////////////////////////////////////////////////////////////
    
    // Create the disk
    for (int p = 0; p < num_disk; p++) {
        struct reb_particle starD = {0};
        starD = reb_tools_orbit_to_particle(r->G, BH, part_mD, part_aD[p], 
            part_eD[p], part_iD[p], part_OD[p], part_oD[p], 
            reb_tools_M_to_f(part_eD[p], part_MD[p]));
	    reb_add(r, starD); 
    }
    /////////////////////////////////////////////////////////////////////////
	// Create the bowl
    if (bowl_type !=0){
        for (int p = 0; p < num_bowl; p++) {
            struct reb_particle Bstar = {0};
            Bstar = reb_tools_orbit_to_particle(r->G, BH, part_mB, part_aB[p], 
                part_eB[p], part_iB[p], part_OB[p], part_oB[p], 
                reb_tools_M_to_f(part_eB[p], part_MB[p]));
            reb_add(r, Bstar); 
        }
    }
	/////////////////////////////////////////////////////////////////////////
    tmax = orbits * orbit_time;
	// Move to center of mass frame
	reb_move_to_com(r);
    /////////////////////////////////////////////////////////////////////////

    if (argc > 1) {
        for (int i=1; i < argc; i++) {
            if ( !strcmp(argv[i], "load_binary") ) {
                i++;
                printf("Loading sim from binary...\n");
                struct reb_simulation* r    = reb_create_simulation_from_binary(argv[i]);
                r->heartbeat                = heartbeat;
                r->visualization            = REB_VISUALIZATION_OPENGL;
                if ( !strcmp(argv[i+1], "num_orbits") && (argc >= 5) ) {
                    i += 2;
                    tmax = r->t + orbit_time*atof(argv[i]);
                }
            }
            else if ( !strcmp(argv[i], "save_loc") ) {
                i++;
                sprintf(tfilename, "%s%s", argv[i], "comp_time.csv");
                ///////////////////////////////////////////////////////////////////
                sprintf(ofilename_a, "%s%s", argv[i], "orbits_disk_a.csv");
                sprintf(ofilename_e, "%s%s", argv[i], "orbits_disk_e.csv");
                sprintf(ofilename_inc, "%s%s", argv[i], "orbits_disk_inc.csv");
                sprintf(ofilename_omega, "%s%s", argv[i], "orbits_disk_omega.csv");
                sprintf(ofilename_Omega, "%s%s", argv[i], "orbits_disk_Omega.csv");
                sprintf(ofilename_M, "%s%s", argv[i], "orbits_disk_M.csv");
                /////////////////////////////////////////////////////////////////////
				sprintf(ofilename1_a, "%s%s", argv[i], "orbits_bowl_a.csv");
                sprintf(ofilename1_e, "%s%s", argv[i], "orbits_bowl_e.csv");
                sprintf(ofilename1_inc, "%s%s", argv[i], "orbits_bowl_inc.csv");
                sprintf(ofilename1_omega, "%s%s", argv[i], "orbits_bowl_omega.csv");
                sprintf(ofilename1_Omega, "%s%s", argv[i], "orbits_bowl_Omega.csv");
                sprintf(ofilename1_M, "%s%s", argv[i], "orbits_bowl_M.csv");
                ////////////////////////////////////////////////////////////////////
                sprintf(efilename, "%s%s", argv[i], "e_and_j.csv");
                sprintf(binprefix, "%s%s", argv[i], "save");
            }   
        }
    }
    else {
        sprintf(tfilename, "%s", "comp_time.csv");
        ///////////////////////////////////////////////////////////////////////////////
        sprintf(ofilename_a, "%s", "orbits_disk_a.csv");
        sprintf(ofilename_e, "%s", "orbits_disk_e.csv");
        sprintf(ofilename_inc, "%s", "orbits_disk_inc.csv");
        sprintf(ofilename_omega, "%s", "orbits_disk_omega.csv");
        sprintf(ofilename_Omega, "%s", "orbits_disk_Omega.csv");
        sprintf(ofilename_M, "%s", "orbits_disk_M.csv");
        //////////////////////////////////////////////////////////////////////////////
		sprintf(ofilename1_a, "%s", "orbits_bowl_a.csv");
        sprintf(ofilename1_e, "%s", "orbits_bowl_e.csv");
        sprintf(ofilename1_inc, "%s", "orbits_bowl_inc.csv");
        sprintf(ofilename1_omega, "%s", "orbits_bowl_omega.csv");
        sprintf(ofilename1_Omega, "%s", "orbits_bowl_Omega.csv");
        sprintf(ofilename1_M, "%s", "orbits_bowl_M.csv");
        /////////////////////////////////////////////////////////////////////////////
        sprintf(efilename, "%s", "e_and_j.csv");
        sprintf(binprefix, "%s", "save");
    }
    // Begin integration
    reb_integrate(r, tmax);	
    /////////////////////////////////////////////////////////////////////////
    printf("\n");
}    

void basic_save(struct reb_simulation* r, char* filename_a, char* filename_e, char* filename_inc, char* filename_omega, char* filename_Omega, char* filename_M, char* filename1_a, char* filename1_e, char* filename1_inc, char* filename1_omega, char* filename1_Omega, char* filename1_M) {
    // disk files
    FILE* of_a = fopen(filename_a, "a");
    FILE* of_e = fopen(filename_e, "a");
    FILE* of_inc = fopen(filename_inc, "a");
    FILE* of_omega = fopen(filename_omega, "a");
    FILE* of_Omega = fopen(filename_Omega, "a");
    FILE* of_M = fopen(filename_M, "a");
    // bowl files
	FILE* of1_a = fopen(filename1_a, "a");
    FILE* of1_e = fopen(filename1_e, "a");
    FILE* of1_inc = fopen(filename1_inc, "a");
    FILE* of1_omega = fopen(filename1_omega, "a");
    FILE* of1_Omega = fopen(filename1_Omega, "a");
    FILE* of1_M = fopen(filename1_M, "a");
    struct reb_particle* BH = reb_get_particle_by_hash(r, reb_hash("BH"));
    for (int p = 0; p < r->N; p++) {
        if (r->particles[p].hash != BH->hash) {
			// save to disk_orbits.csv
			if (p <= num_disk) {	
				struct reb_orbit part_orbit = reb_tools_particle_to_orbit(r->G, r->particles[p], *BH);
				fprintf(of_a, "%e", part_orbit.a);
				fprintf(of_a, ",");
                fprintf(of_e, "%e", part_orbit.e);
				fprintf(of_e, ",");
                fprintf(of_inc, "%e", part_orbit.inc);
				fprintf(of_inc, ",");
                fprintf(of_omega, "%e", part_orbit.omega);
				fprintf(of_omega, ",");
                fprintf(of_Omega, "%e", part_orbit.Omega);
				fprintf(of_Omega, ",");
                fprintf(of_M, "%e", part_orbit.M);
				fprintf(of_M, ",");
			}
			// save to bowl_orbits.csv
			else if((p > num_disk) && (bowl_type != 0)){	
				struct reb_orbit part_orbit = reb_tools_particle_to_orbit(r->G, r->particles[p], *BH);
                fprintf(of1_a, "%e", part_orbit.a);
				fprintf(of1_a, ",");
                fprintf(of1_e, "%e", part_orbit.e);
				fprintf(of1_e, ",");
                fprintf(of1_inc, "%e", part_orbit.inc);
				fprintf(of1_inc, ",");
                fprintf(of1_omega, "%e", part_orbit.omega);
				fprintf(of1_omega, ",");
                fprintf(of1_Omega, "%e", part_orbit.Omega);
				fprintf(of1_Omega, ",");
                fprintf(of1_M, "%e", part_orbit.M);
				fprintf(of1_M, ",");
				}
        	}
    	}
    fprintf(of_a, ",%e\n", r->t/(2.*M_PI));
    fprintf(of1_a, ",%e\n", r->t/(2.*M_PI));
    fprintf(of_e, ",%e\n", r->t/(2.*M_PI));
    fprintf(of1_e, ",%e\n", r->t/(2.*M_PI));
    fprintf(of_inc, ",%e\n", r->t/(2.*M_PI));
    fprintf(of1_inc, ",%e\n", r->t/(2.*M_PI));
    fprintf(of_omega, ",%e\n", r->t/(2.*M_PI));
    fprintf(of1_omega, ",%e\n", r->t/(2.*M_PI));
    fprintf(of_Omega, ",%e\n", r->t/(2.*M_PI));
    fprintf(of1_Omega, ",%e\n", r->t/(2.*M_PI));
    fprintf(of_M, ",%e\n", r->t/(2.*M_PI));
    fprintf(of1_M, ",%e\n", r->t/(2.*M_PI));
	fclose(of_a);
    fclose(of1_a);
    fclose(of_e);
    fclose(of1_e);
    fclose(of_inc);
    fclose(of1_inc);
    fclose(of_omega);
    fclose(of1_omega);
    fclose(of_Omega);
    fclose(of1_Omega);
    fclose(of_M);
    fclose(of1_M);

}

void heartbeat(struct reb_simulation* r) {
    double last_dt = r->dt_last_done;
    orbit_counter += last_dt;
    save_counter += last_dt;
    if ((orbit_counter >= orbit_time/num_saves_per_orbit) || (orbit_counter == 0.)) {
        // Output progress to screen
        reb_output_timing(r, tmax);
        /////////////////////////////////////////////////////////////////////
        // Record the time between heartbeats in milliseconds
        int msec_since_last_hb = (clock() - last_hb)*1000/CLOCKS_PER_SEC;
        last_hb = clock();
        FILE* tf = fopen(tfilename, "a");
        fprintf(tf, "%d,%e\n", msec_since_last_hb, r->t/(2.*M_PI));
        fclose(tf);
        /////////////////////////////////////////////////////////////////////
        // Save orbits to a file
        //clock_t start = clock(), diff;
		printf("Save orbits...\n");
        basic_save(r, ofilename_a, ofilename_e, ofilename_inc, ofilename_omega, ofilename_Omega, ofilename_M, ofilename1_a, ofilename1_e, ofilename1_inc, ofilename1_omega, ofilename1_Omega, ofilename1_M);
        printf("Done saving...\n");
        /////////////////////////////////////////////////////////////////////
        // Save energy and ang momentum to file
        printf("Save e_and_j...\n");
        FILE* ef = fopen(efilename, "a");
        struct reb_vec3d j = reb_tools_angular_momentum(r);
        double total_J = sqrt(j.x*j.x + j.y*j.y + j.z*j.z);
        fprintf(ef, "%.10e,%.10e,%e\n", reb_tools_energy(r), total_J, r->t/(2.*M_PI));
        fclose(ef);
        printf("Done saving...\n");
        /////////////////////////////////////////////////////////////////////
        // Reset counter
        orbit_counter -= orbit_time/num_saves_per_orbit;
    if ((save_counter >= tmax/10) || (save_counter == 0.)) {
        // Save binary of the simulation
        char bin_filename[100];
        sprintf(bin_filename, "%s_%d.bin", binprefix, bin_save_count);
        reb_output_binary(r, bin_filename);
        bin_save_count++;
        /////////////////////////////////////////////////////////////////////
        // Reset counter
        save_counter = 0.;
    	}
	}
}