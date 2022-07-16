#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 18:54:23 2022

@author: bayat2
"""

        """for enum_iter in np.arange(self.I["num_ietartion"]):
            done_flag=False
            self.I["num_seg_iters"][enum_iter]=round(self.I["num_time_p_iters"][enum_iter]/20)+1
            if enum_iter>1:
                if self.I["num_time_p_iters"][enum_iter]==self.I["num_time_p_iters"][enum_iter-1]:
                    self.I["num_nodes_each_seg_iters"][enum_iter]=self.I["num_nodes_each_seg_iters"][enum_iter-1]
                    done_flag=True
            
            if done_flag==False:
                num_LGR_nodes=round(self.I["num_time_p_iters"][enum_iter]/self.I["num_seg_iters"][enum_iter])
                last_node=LGR_Class(num_LGR_nodes,False).LGR_Nodes[-1] 
                while last_node>1-2/(self.I["num_time_p_iters"][enum_iter]/self.I["num_seg_iters"][enum_iter]):
                    num_LGR_nodes=num_LGR_nodes-1
                    last_node=LGR_Class(num_LGR_nodes,False).LGR_Nodes[-1] 
                self.I["num_nodes_each_seg_iters"][enum_iter]=num_LGR_nodes"""