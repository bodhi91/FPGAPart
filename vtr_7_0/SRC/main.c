/**
 VPR is a CAD tool used to conduct FPGA architecture exploration.  It takes, as input, a technology-mapped netlist and a description of the FPGA architecture being investigated.  
 VPR then generates a packed, placed, and routed FPGA (in .net, .place, and .route files respectively) that implements the input netlist.
 
 This file is where VPR starts execution.

 Key files in VPR:
 1.  libarchfpga/physical_types.h - Data structures that define the properties of the FPGA architecture
 2.  vpr_types.h - Very major file that defines the core data structures used in VPR.  This includes detailed architecture information, user netlist data structures, and data structures that describe the mapping between those two.
 3.  globals.h - Defines the global variables used by VPR.
 */

#include <cstdio>
#include <cstring>
#include <ctime>
using namespace std;

#include "vpr_api.h"

// operates on logical blocks, which are in the BLIF input to VPR

static void print_triton_part_logical_hypergraph(bool hcrt = false, bool crt = false, bool nname = false) {
  FILE *hypergraph;
  hypergraph = fopen("logical_hypergraph_tpart.txt", "w");
  int wt_status;
  if (hcrt == true && crt == false) {
	wt_status = 1;
  } else if (hcrt == false && crt == true) {
	wt_status = 10;
  } else if (hcrt == true && crt == true) {
	wt_status = 11;
  } else {
	wt_status = -1;
  }	

  fprintf(hypergraph, "%d %d", num_logical_nets, num_logical_blocks);
  if (wt_status != 1) {
	fprintf(hypergraph, " %d\n", wt_status);
  } else {
	fprintf(hypergraph, "\n");
  }
  for (int i = 0; i < num_logical_nets; i++) {
	if (hcrt == true) {
		fprintf(hypergraph, "%f ", logical_net_criticalities[i]);
	}
	for(int j = 0; j <= vpack_net[i].num_sinks; j++) {
	  fprintf(hypergraph, "%d ", vpack_net[i].node_block[j]+1);
	}
	if (nname == true) {
		fprintf(hypergraph, "//%s", vpack_net[i].name);
	}
	fprintf(hypergraph, "\n");
  }
  if (crt == true) {
	for (int i = 0; i < num_logical_blocks; i++) {
	  fprintf(hypergraph, "%f ", logical_block_criticalities[i]);
	  for (int j = 0; j < logical_block_types[i].size(); j++) {
		fprintf(hypergraph, "%d ", logical_block_types[i][j]);
	  }
	  fprintf(hypergraph, "\n");
	}
  } else {
	for (int i = 0; i < num_logical_blocks; i++) {
	  for (int j = 0; j < logical_block_types[i].size(); j++) {
		fprintf(hypergraph, "%d ", logical_block_types[i][j]);
	  }
	  fprintf(hypergraph, "\n");
	}
  }
  fclose(hypergraph);

  FILE* blocks_file = fopen("logical_blocks.txt", "w");
  for(int i = 0; i < num_logical_blocks; i++){
	fprintf(blocks_file, "%s %s\n", logical_block[i].name, logical_block[i].model->name);
  }
  fclose(blocks_file);	
}

static void print_hmetis_logical_graph(bool hwt = false, bool vwt = false, bool nname = false) {
  FILE *hypergraph;
  hypergraph = fopen("logical_hypergraph.txt", "w");
  int wt_status;
  if (hwt == true && vwt == false) {
	wt_status = 1;
  } else if (hwt == false && vwt == true) {
	wt_status = 10;
  } else if (hwt == true && vwt == true) {
	wt_status = 11;
  } else {
	wt_status = -1;
  }	

  fprintf(hypergraph, "%d %d", num_logical_nets, num_logical_blocks);
  if (wt_status != 1) {
	fprintf(hypergraph, " %d\n", wt_status);
  } else {
	fprintf(hypergraph, "\n");
  }
  for (int i = 0; i < num_logical_nets; i++) {
	if (hwt == true) {
		fprintf(hypergraph, "%f ", logical_net_criticalities[i]);
	}
    for(int j = 0; j <= vpack_net[i].num_sinks; j++) {
      fprintf(hypergraph, "%d ", vpack_net[i].node_block[j]+1);
    }
    if (nname == true) {
		fprintf(hypergraph, "//%s", vpack_net[i].name);
	}
    fprintf(hypergraph, "\n");
  }
  if (vwt == true) {
	for (int i = 0; i < num_logical_blocks; i++) {
	  fprintf(hypergraph, "%f\n", logical_block_criticalities[i]);
	}
  }
  fclose(hypergraph);

  FILE* blocks_file = fopen("logical_blocks.txt", "w");
  for(int i = 0; i < num_logical_blocks; i++){
    fprintf(blocks_file, "%s %s\n", logical_block[i].name, logical_block[i].model->name);
  }
  fclose(blocks_file);

}

// operates on CLBs, which are the complex blocks that are placed during placement
static void print_hmetis_clb_graph() {
			/* ANDRE: printing the hypergraph */
			FILE *hypergraph;
			hypergraph = fopen("clb_hypergraph.txt", "w");

			fprintf(hypergraph, "%d %d\n", num_nets, num_blocks); // first line of file: numnets, numblocks, format (weighted edges and nodes)
			for(int i = 0; i < num_nets; i++){
			  //fprintf(hypergraph, "%d", 1 - clb_net[i].is_global); //weight of the edge, 0 if it is global
			  for(int j = 0; j <= clb_net[i].num_sinks; j++)
			    fprintf(hypergraph, "%d ", clb_net[i].node_block[j]+1);
			  fprintf(hypergraph, "//%s", clb_net[i].name);
			  fprintf(hypergraph, "\n");
			}
			fclose(hypergraph);

			FILE* blocks_file = fopen("clb_blocks.txt", "w");
			for(int i = 0; i < num_blocks; i++){
			  fprintf(blocks_file, "%s %s\n", block[i].name, block[i].pb->pb_graph_node->pb_type->name);
			}
			fclose(blocks_file);


}

/**
 * VPR program
 * Generate FPGA architecture given architecture description
 * Pack, place, and route circuit into FPGA architecture
 * Electrical timing analysis on results
 *
 * Overall steps
 * 1.  Initialization
 * 2.  Pack
 * 3.  Place-and-route and timing analysis
 * 4.  Clean up
 */
int main(int argc, char **argv) {
	t_options Options;
	t_arch Arch;
	t_vpr_setup vpr_setup;
	clock_t entire_flow_begin,entire_flow_end;

	entire_flow_begin = clock();
	try{
		/* Read options, architecture, and circuit netlist */
		vpr_init(argc, argv, &Options, &vpr_setup, &Arch);

		/* If the user requests packing, do packing */
		if (vpr_setup.PackerOpts.doPacking) {
			vpr_pack(vpr_setup, Arch);
			vpr_init_pre_place_and_route(vpr_setup, Arch);
         	        print_hmetis_logical_graph(true, true);
					print_triton_part_logical_hypergraph(true, true);
         	        print_hmetis_clb_graph();
		} else {
			vpr_init_pre_place_and_route(vpr_setup, Arch);
                }

		if (vpr_setup.PlacerOpts.doPlacement || vpr_setup.RouterOpts.doRouting) {
			
			#ifdef INTERPOSER_BASED_ARCHITECTURE
			vpr_setup_interposer_cut_locations(Arch);
			#endif
			
			vpr_place_and_route(vpr_setup, Arch);
	#if 0
			if(vpr_setup.RouterOpts.doRouting) {
				vpr_resync_post_route_netlist_to_TI_CLAY_v1_architecture(&Arch);
			}
	#endif
		}

		if (vpr_setup.PowerOpts.do_power) {
			vpr_power_estimation(vpr_setup, Arch);
		}
	
		entire_flow_end = clock();
	
		#ifdef CLOCKS_PER_SEC
			vpr_printf_info("The entire flow of VPR took %g seconds.\n", 
					(float)(entire_flow_end - entire_flow_begin) / CLOCKS_PER_SEC);
		#else
			vpr_printf_info("The entire flow of VPR took %g seconds.\n", 
					(float)(entire_flow_end - entire_flow_begin) / CLK_PER_SEC);
		#endif
	
		/* free data structures */
		vpr_free_all(Arch, Options, vpr_setup);
	}
	catch(t_vpr_error* vpr_error){
		vpr_print_error(vpr_error);
	}
	
	/* Return 0 to single success to scripts */
	return 0;
}




