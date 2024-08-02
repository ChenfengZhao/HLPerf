import simpy
from utils import *
import timeit
import argparse

def read_nod_src(env, fifo_out_nod_src_list, nidx_begin, nidx_end, L, II, L_mem, deg_list):
    """read nod_src from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_out_nod_src_list : list of simpy.store
        The list of output FIFO of neighbor index range (or degree) of each node . length is 6
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    deg_list : list of int

    """
    print("Simulation starts at %d" % env.now)

    for n in range(nidx_begin, nidx_end):
        # yield env.timeout(L_mem) # memory latency
        # yield env.timeout(15) # memory latency

        if n == nidx_begin:
            yield env.timeout(L+L_mem)
        else:
            yield env.timeout(II)

        for fifo_idx in range(len(fifo_out_nod_src_list)):
            yield fifo_out_nod_src_list[fifo_idx].put(deg_list[n])

    print("read_nod_src ends at %d" % env.now)

def read_edge_src(env, fifo_in_nod_src, fifo_out_tmp_src, nidx_begin, nidx_end, L, II, L_mem):
    """read edge_src from memory
    
    Instead of putting actual edge index into the FIFO, we put fake values

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_tmp_src : simpy.store
        Output FIFO of neibhbor indices of each node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        # yield env.timeout(L_mem) # memory latency
        # yield env.timeout(15) # memory latency

        for e in range(deg):
            if e == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_tmp_src.put(1)

    print("read_edge_src ends at %d" % env.now)

def read_feat_in_agg(env, fifo_in_tmp_src, fifo_in_nod_src, fifo_out_ft_in_agg, nidx_begin, nidx_end, L, II, L_mem):
    """read input features to be aggregated from memory

    Instead of putting acutal input feature of nbrs, we put feke values

    Parameters  
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_tmp_src :simpy.store
        Input FIFO of neighbor indices of each node
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_ft_in_agg : simpy.store
        Output FIFO of input feature to be aggretated
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for e in range(deg):
            yield env.timeout(3)
            _ = yield fifo_in_tmp_src.get()

            # yield env.timeout(L_mem) # memory latency
            # yield env.timeout(15) # memory latency

            for i in range(2):
                if i == 0:
                    yield env.timeout(L)
                else:
                    yield env.timeout(II)

                yield fifo_out_ft_in_agg.put(1)
                
    
    print("read_feat_in_agg ends at %d" % env.now)

def agg_feat_in(env, fifo_in_ft_in_agg, fifo_in_nod_src, fifo_out_ft_h_agg, nidx_begin, nidx_end, L_list, II_list):
    """aggregate features

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_in_agg : simpy.store
        Input FIFO of input feature to be aggretated
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_ft_h_agg : simpy.store
        Output FIFO of aggretaed fature
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency of the loop
    II_list : list of int
        The II of the loop
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        if deg != 0:
            for e in range(deg):

                for i in range(2):
                    if e == 0 and i == 0:
                        yield env.timeout(L_list[0])
                    else:
                        yield env.timeout(II_list[0])
                    _ = yield fifo_in_ft_in_agg.get()

            for i in range(2):
                if i == 0:
                    yield env.timeout(L_list[1])
                else:
                    yield env.timeout(II_list[1])
            
                yield fifo_out_ft_h_agg.put(1)
    
    print("agg_feat_in ends at %d" % env.now)

def read_feat_in_tar(env, fifo_out_ft_in_tar, nidx_begin, nidx_end, L, II, L_mem):
    """read the target feature from memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_out_ft_in_tar : simpy.store
        Output fifo of input feature of the target node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory read access
    """

    for n in range(nidx_begin, nidx_end):
        # yield env.timeout(L_mem)

        for i in range(2):
            if n == nidx_begin and i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            yield fifo_out_ft_in_tar.put(1)

    print("read_feat_in_tar ends at %d" % env.now)

def concat_rst(env, fifo_in_ft_h_agg, fifo_in_ft_in_tar, fifo_in_nod_src, fifo_out_rst_cat, nidx_begin, nidx_end, L, II):
    """combine/concat the aggregated results ft_h_agg and target features ft_in_tar_stream: (1 + eps)*ft_in_tar_stream[i] + ft_h_agg[i]

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_ft_h_agg : simpy.store
        Input FIFO of result aggregated feature
    fifo_in_ft_in_tar : simpy.store
        Input FIFO of result target feature
    fifo_in_nod_src : simpy.store
        Input FIFO of neighbor index range (or degree) of each node 
    fifo_out_rst_cat : simpy.store
        Ouput FIFO of result cat feature
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    """

    for n in range(nidx_begin, nidx_end):
        deg = yield fifo_in_nod_src.get()
        yield env.timeout(3)

        for i in range(2):
            if i == 0:
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            if deg != 0:
                _ = yield fifo_in_ft_in_tar.get()
                _ = yield fifo_in_ft_h_agg.get()
                yield fifo_out_rst_cat.put(1)
            else:
                _ = yield fifo_in_ft_in_tar.get()
                yield fifo_out_rst_cat.put(1)
    
    print("concat_rst ends at %d" % env.now)

def mlp0(env, fifo_in_rst_cat, fifo_out_rst_mlp0_p1, nidx_begin, nidx_end, L_list, II_list):
    """MLP layer 0 update/apply phase for the target feature (phase 1)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_cat : simpy.store
        Input FIFO of result cat feature
    fifo_out_rst_mlp0_p1 : simpy.store
        Output FIFO of result feature (phase 1) of each target node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    """

    for n in range(nidx_begin, nidx_end):

        for i in range(2):
            if i == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])
            _ = yield fifo_in_rst_cat.get()
        
        yield env.timeout(L_list[1] + (32-1) * II_list[1])

        for kd in range(8):
            if kd == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
            
            yield fifo_out_rst_mlp0_p1.put(1)

    print("mlp0 ends at %d" % env.now)

def mlp0_sum(env, fifo_in_rst_mlp0_p1, fifo_out_rst_mlp0, nidx_begin, nidx_end, L_list, II_list):
    """MLP layer 0 update/apply phase for the target feature (phase 2)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_mlp0_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each target node
    fifo_out_rst_mlp0 : simpy.store
        Output FIFO of result feature (phase 2) of each target node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    """

    for n in range(nidx_begin, nidx_end):
        for kd in range(8):
            if kd == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])
            _ = yield fifo_in_rst_mlp0_p1.get()
        
        yield env.timeout(L_list[1] + (32-1)*II_list[1])

        for i in range(2):
            if i == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
            yield fifo_out_rst_mlp0.put(1)
    
    print("mlp0_sum ends at %d" % env.now)

def mlp1(env, fifo_in_rst_mlp0, fifo_out_rst_mlp1_p1, nidx_begin, nidx_end, L_list, II_list):
    """MLP layer 1 update/apply phase for the target feature (phase 1)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_mlp0 : simpy.store
        Input FIFO of MLP0 result feature (phase 2) of each target node
    fifo_out_rst_mlp1_p1 : simpy.store
        Output FIFO of result feature (phase 1) of each target node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    """

    for n in range(nidx_begin, nidx_end):

        for i in range(2):
            if i == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])
            _ = yield fifo_in_rst_mlp0.get()
        
        yield env.timeout(L_list[1] + (32-1) * II_list[1])

        for kd in range(8):
            if kd == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
            
            yield fifo_out_rst_mlp1_p1.put(1)

    print("mlp1 ends at %d" % env.now)

def mlp1_sum(env, fifo_in_rst_mlp1_p1, fifo_out_rst_mlp1, nidx_begin, nidx_end, L_list, II_list):
    """MLP layer 1 update/apply phase for the target feature (phase 2)

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_mlp1_p1 : simpy.store
        Input FIFO of result feature (phase 1) of each target node
    fifo_out_rst_mlp1 : simpy.store
        Output FIFO of result feature (phase 2) of each target node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L_list : list of int
        The latency list of loops
    II_list : list of int
        The II list of loops
    """

    for n in range(nidx_begin, nidx_end):
        for kd in range(8):
            if kd == 0:
                yield env.timeout(L_list[0])
            else:
                yield env.timeout(II_list[0])
            _ = yield fifo_in_rst_mlp1_p1.get()
        
        yield env.timeout(L_list[1] + (32-1)*II_list[1])

        for i in range(2):
            if i == 0:
                yield env.timeout(L_list[2])
            else:
                yield env.timeout(II_list[2])
            yield fifo_out_rst_mlp1.put(1)
    
    print("mlp1_sum ends at %d" % env.now)

def write_mlp1_rst_mem(env, fifo_in_rst_mlp1, nidx_begin, nidx_end, L, II, L_mem):
    """write results of MLP layer 1 to memory

    Parameters
    ----------
    env : simpy.env
        The enviornment of simulation
    fifo_in_rst_mlp1 : simpy.store
        Input FIFO of MLP 1 result feature of each target node
    nidx_begin : int
        The begin node index
    nidx_end : int
        The end node index
    L : int
        The latency of the loop
    II : int
        The II of the loop
    L_mem : int
        the latency of memory write access
    """

    for n in range(nidx_begin, nidx_end):

        for i in range(2):
            if i == 0:
                # yield env.timeout(L+L_mem)
                yield env.timeout(L)
            else:
                yield env.timeout(II)
            
            _ = yield fifo_in_rst_mlp1.get()

    print("write_mlp1_rst_mem ends at %d" % env.now)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_option(parser)
    args = parser.parse_args()

    # read input graph
    format = "gnnhls"
    # input_fn = "./test/csr_indptr_trans.txt"
    # input_fn = "./test/csr_indptr_trans_v2.txt"

    # input_fn = "./data/OGBG-MOLTOX21/csr_indptr_trans.txt"
    # input_fn = "./data/OGBG-MOLHIV/csr_indptr_trans.txt"
    # input_fn = "./data/OGBN-ARXIV/csr_indptr_trans.txt"
    # input_fn = "./data/OGBN-PROTEINS/csr_indptr_trans.txt"

    # input_fn = input_fn_dict[args.ds_fn]
    input_fn = gen_input_fn(args.ds_fn)
    print("Input dataset:", args.ds_fn, ":", input_fn)

    # hardware parameters (MHz)
    freq = 298.9
    print("CLK Frequency:", freq, "MHz")

    # generate degree list of the input graph
    deg_list = load_graph_deg(input_fn, format)
    # deg_list = [2] * 5000 + [50] * 5000
    # deg_list = [2, 50] * 5000

    # calculate node number and edge number
    node_num = len(deg_list)
    edge_num = sum(deg_list)
    print("node #:", node_num, "edge #:", edge_num)

    nidx_begin = 0
    nidx_end = node_num

    start = timeit.default_timer()

    env = simpy.Environment()

    # read nod_src from memory
    fifo_nod_src_list = [simpy.Store(env, capacity=4), simpy.Store(env, capacity=5), simpy.Store(env, capacity=6), simpy.Store(env, capacity=8)]
    L = 4
    II = 1
    L_mem = 64
    env.process(read_nod_src(env, fifo_nod_src_list, nidx_begin, nidx_end, L, II, L_mem, deg_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = L_mem + L + II * (node_num - 1)
    print("Estimated exec time of read_nod_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read edge_src from memory
    fifo_tmp_src = simpy.Store(env, capacity=10)
    L = 75
    II = 1
    L_mem = 64
    env.process(read_edge_src(env, fifo_nod_src_list[0], fifo_tmp_src, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3) * node_num + L * node_num + II * (edge_num - node_num)
    print("Estimated exec time of read_edge_src alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # read input features to be aggregated from memory
    fifo_ft_in_agg = simpy.Store(env, capacity=10)
    L = 74
    II = 1
    L_mem = 64
    env.process(read_feat_in_agg(env, fifo_tmp_src, fifo_nod_src_list[1], fifo_ft_in_agg, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + (1) * edge_num + (L + II) * edge_num
    print("Estimated exec time of read_feat_in_agg alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # aggregate features
    fifo_ft_h_agg = simpy.Store(env, capacity=10)
    L_list = [8, 2]
    II_list = [3, 1]
    env.process(agg_feat_in(env, fifo_ft_in_agg, fifo_nod_src_list[2], fifo_ft_h_agg, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = 3 * node_num + L_list[0] * node_num + II_list[0] * edge_num + (1 + L_list[1] + II_list[1]) * node_num
    print("Estimated exec time of agg_feat_in alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))
    
    # read the target feature from memory
    fifo_ft_in_tar = simpy.Store(env, capacity=10)
    L = 76
    II = 1
    L_mem = 64
    env.process(read_feat_in_tar(env, fifo_ft_in_tar, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_mem + L + II) * node_num
    print("Estimated exec time of read_feat_in_tar alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # combine/concat the aggregated results ft_h_agg and target features ft_in_tar_stream: (1 + eps)*ft_in_tar_stream[i] + ft_h_agg[i]
    fifo_rst_cat = simpy.Store(env, capacity=10)
    L = 14
    II = 1
    env.process(concat_rst(env, fifo_ft_h_agg, fifo_ft_in_tar, fifo_nod_src_list[3], fifo_rst_cat, nidx_begin, nidx_end, L, II))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (3 + L + II) * node_num
    print("Estimated exec time of concat_rst alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # MLP layer 0 update/apply phase for the target feature (phase 1)
    fifo_rst_mlp0_p1 = simpy.Store(env, capacity=16)
    L_list = [2, 14, 2]
    II_list = [1, 2, 1]
    env.process(mlp0(env, fifo_rst_cat, fifo_rst_mlp0_p1, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + II_list[0] + L_list[1] + (32-1) * II_list[1] + L_list[2] + (8 - 1) * II_list[2]) * node_num
    print("Estimated exec time of mlp0 alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # MLP layer 0 update/apply phase for the target feature (phase 2)
    fifo_rst_mlp0 = simpy.Store(env, capacity=10)
    L_list = [2, 9, 4]
    II_list = [1, 2, 1]
    env.process(mlp0_sum(env, fifo_rst_mlp0_p1, fifo_rst_mlp0, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + (8 - 1) * II_list[0] + L_list[1] + (32-1)*II_list[1] + L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of mlp0_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # MLP layer 1 update/apply phase for the target feature (phase 1)
    fifo_rst_mlp1_p1 = simpy.Store(env, capacity=16)
    L_list = [2, 14, 2]
    II_list = [1, 2, 1]
    env.process(mlp1(env, fifo_rst_mlp0, fifo_rst_mlp1_p1, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + II_list[0] + L_list[1] + (32-1) * II_list[1] + L_list[2] + (8 - 1) * II_list[2]) * node_num
    print("Estimated exec time of mlp1 alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # MLP layer 1 update/apply phase for the target feature (phase 2)
    fifo_rst_mlp1 = simpy.Store(env, capacity=10)
    L_list = [2, 9, 4]
    II_list = [1, 2, 1]
    env.process(mlp1_sum(env, fifo_rst_mlp1_p1, fifo_rst_mlp1, nidx_begin, nidx_end, L_list, II_list))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L_list[0] + (8 - 1) * II_list[0] + L_list[1] + (32-1)*II_list[1] + L_list[2] + II_list[2]) * node_num
    print("Estimated exec time of mlp1_sum alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))

    # write results of MLP layer 1 to memory
    L = 73
    II = 1
    L_mem = 64
    env.process(write_mlp1_rst_mem(env, fifo_rst_mlp1, nidx_begin, nidx_end, L, II, L_mem))
    # Estimated ideal execution time of each function alone, without other functions and stalls
    exec_time_alone = (L + II) * node_num
    print("Estimated exec time of write_mlp1_rst_mem alone: %d cycles," % exec_time_alone, "%f s" % (exec_time_alone * (1000/freq) * 1e-9))


    # run the simulation
    env.run()

    # measure the actural running time (seconds)
    end = timeit.default_timer()
    print("Actual running time: %s Seconds"%(end-start))

    # measure the simulation resultes (cycles)
    print("Simulation ends at %d" % env.now)
    print("Calculated simulation time (s): ", (env.now * 1000 / freq * 1e-9))

    
