f_2_timer.start();
    for (int i=0; i<N; i++) {
         //--  onsites(EVEN) {
        //--              f_2[X] = f_2[X + e_z]*f_2[X - e_z];
        //--          }
        {
          Field<SU<2, double>> & HILA_field_f_2 = f_2;
          HILA_field_f_2.check_alloc();
          if (HILA_field_f_2.fs->mylattice.ptr() != lattice.ptr()) {hila::out0 << "File /home/haaaaron/HILA/applications/matrix-test/src/matrix-test.cpp on line 36: Field f_2 initialized on different lattice\n";
            hila::terminate(1);
          }
          if (!HILA_field_f_2.is_initialized(opp_parity(Parity::even))) {
            hila::out0 << "File /home/haaaaron/HILA/applications/matrix-test/src/matrix-test.cpp on line 36:\n Variable f_2 is not properly initialized\n";
            hila::terminate(1);
          }
          dir_mask_t  _dir_mask_ = 0;
          auto [mask1, send_params_1] = HILA_field_f_2.start_gather_test(e_z, Parity::even);
          auto [mask2, send_params_2] = HILA_field_f_2.start_gather_test(- e_z, Parity::even);
          _dir_mask_ |= mask1;
          _dir_mask_ |= mask2;
          HILA_field_f_2.send_buffers(send_params_1);
          HILA_field_f_2.send_buffers(send_params_2);
          const lattice_struct & loop_lattice = lattice.ref();
          const int loop_begin = loop_lattice.loop_begin(Parity::even);
          const int loop_end   = loop_lattice.loop_end(Parity::even);
          for (int _wait_i_ = 0; _wait_i_ < 2; ++_wait_i_) {
            for(int HILA_index = loop_begin; HILA_index < loop_end; ++HILA_index) {
              if (((loop_lattice.wait_arr_[HILA_index] & _dir_mask_) != 0) == _wait_i_) {
                const SU<2, double> HILA_field_f_2_dir1 = HILA_field_f_2.get_value_at(HILA_field_f_2.fs->neighbours[e_z][HILA_index]);
                const SU<2, double> HILA_field_f_2_dir2 = HILA_field_f_2.get_value_at(HILA_field_f_2.fs->neighbours[- e_z][HILA_index]);
                SU<2, double> HILA_field_f_2_at_X;
                // Initial value of variable HILA_field_f_2_at_X not needed
                {
                  HILA_field_f_2_at_X = HILA_field_f_2_dir1*HILA_field_f_2_dir2;
                }
                HILA_field_f_2.set_value_at(HILA_field_f_2_at_X, HILA_index);
              }
            }
            if (_dir_mask_ == 0) break;    // No need for another round
            HILA_field_f_2.wait_gather(e_z, Parity::even);
            HILA_field_f_2.wait_gather(- e_z, Parity::even);
          }
          hila::set_allreduce(true);
          HILA_field_f_2.mark_changed(Parity::even);
        }
        //----------
        
         //--  onsites(ODD) {
        //--              f_2[X] = f_2[X + e_z]*f_2[X - e_z];
        //--          }
        {
          Field<SU<2, double>> & HILA_field_f_2 = f_2;
          HILA_field_f_2.check_alloc();
          if (HILA_field_f_2.fs->mylattice.ptr() != lattice.ptr()) {hila::out0 << "File /home/haaaaron/HILA/applications/matrix-test/src/matrix-test.cpp on line 39: Field f_2 initialized on different lattice\n";
            hila::terminate(1);
          }
          if (!HILA_field_f_2.is_initialized(opp_parity(Parity::odd))) {
            hila::out0 << "File /home/haaaaron/HILA/applications/matrix-test/src/matrix-test.cpp on line 39:\n Variable f_2 is not properly initialized\n";
            hila::terminate(1);
          }
          dir_mask_t  _dir_mask_ = 0;
          auto [mask1, send_params_1] = HILA_field_f_2.start_gather_test(e_z, Parity::odd);
          auto [mask2, send_params_2] = HILA_field_f_2.start_gather_test(- e_z, Parity::odd);
          _dir_mask_ |= mask1;
          _dir_mask_ |= mask2;
          HILA_field_f_2.send_buffers(send_params_1);
          HILA_field_f_2.send_buffers(send_params_2);
          const lattice_struct & loop_lattice = lattice.ref();
          const int loop_begin = loop_lattice.loop_begin(Parity::odd);
          const int loop_end   = loop_lattice.loop_end(Parity::odd);
          for (int _wait_i_ = 0; _wait_i_ < 2; ++_wait_i_) {
            for(int HILA_index = loop_begin; HILA_index < loop_end; ++HILA_index) {
              if (((loop_lattice.wait_arr_[HILA_index] & _dir_mask_) != 0) == _wait_i_) {
                const SU<2, double> HILA_field_f_2_dir1 = HILA_field_f_2.get_value_at(HILA_field_f_2.fs->neighbours[e_z][HILA_index]);
                const SU<2, double> HILA_field_f_2_dir2 = HILA_field_f_2.get_value_at(HILA_field_f_2.fs->neighbours[- e_z][HILA_index]);
                SU<2, double> HILA_field_f_2_at_X;
                // Initial value of variable HILA_field_f_2_at_X not needed
                {
                  HILA_field_f_2_at_X = HILA_field_f_2_dir1*HILA_field_f_2_dir2;
                }
                HILA_field_f_2.set_value_at(HILA_field_f_2_at_X, HILA_index);
              }
            }
            if (_dir_mask_ == 0) break;    // No need for another round
            HILA_field_f_2.wait_gather(e_z, Parity::odd);
            HILA_field_f_2.wait_gather(- e_z, Parity::odd);
          }
          hila::set_allreduce(true);
          HILA_field_f_2.mark_changed(Parity::odd);
        }
        //----------
        
    }
    hila::synchronize();
    f_2_timer.stop();