    init_global_streams();
    gpuEvent_t halo_event_2;
    gpuEventCreate(&halo_event_2);
    f_2_timer.start();
    roctxRangePush("Loop Start");
    for (int i=0; i<1000; i++) {
             //--  onsites(EVEN) {
            //--                  for (int i=0; i<N; i++)
            //--                          f_2[X] = f_2[X]*f_2[X + e_z]*f_2[X - e_z].dagger();
            //--              }
            {
              Field<SU<3, float>> & HILA_field_f_2 = f_2;
              HILA_field_f_2.check_alloc();
              if (HILA_field_f_2.fs->mylattice.ptr() != lattice.ptr()) {hila::out0 << "File /projappl/project_462000949/lumi-po-benchmark/benchmarks/matrix-test/src/matrix-test.cpp on line 35: Field f_2 initialized on different lattice\n";
                hila::terminate(1);
              }
              if (!HILA_field_f_2.is_initialized(ALL)) {
                hila::out0 << "File /projappl/project_462000949/lumi-po-benchmark/benchmarks/matrix-test/src/matrix-test.cpp on line 35:\n Variable f_2 is not properly initialized\n";
                hila::terminate(1);
              }
              dir_mask_t  _dir_mask_ = 0;
              roctxRangePush("Start pack");
              
              auto [var_1, mpi_1]= HILA_field_f_2.start_gather_test(e_z, Parity::even);
              auto [var_2, mpi_2]= HILA_field_f_2.start_gather_test(- e_z, Parity::even);
              gpuEventRecord(halo_event_2, halo_stream);
              roctxRangePop();
              roctxRangePush("Kernel bulk");
              const lattice_struct & loop_lattice = lattice.ref();
              int _hila_loop_begin = loop_lattice.loop_begin(Parity::even);
              int _hila_loop_end = loop_lattice.loop_end(Parity::even);
              int N_blocks = (_hila_loop_end - _hila_loop_begin + N_threads - 1)/N_threads;
              hipLaunchKernelGGL(HILA_kernel_main_K31, dim3(N_blocks), dim3(N_threads), 0, bulk_stream, _hila_loop_begin, _hila_loop_end, HILA_field_f_2.fs->payload);

              roctxRangePop();
              roctxRangePush("MPI send");
              gpuStreamWaitEvent(halo_stream, halo_event_2, 0);
              perform_mpi_send(mpi_1);
              perform_mpi_send(mpi_2);

              roctxRangePop();
              roctxRangePush("Wait and unpack");
              HILA_field_f_2.wait_gather(e_z, Parity::even);
              HILA_field_f_2.wait_gather(- e_z, Parity::even);
              roctxRangePop();
              check_device_error("HILA_kernel_main_K31");
              hila::set_allreduce(true);
              HILA_field_f_2.mark_changed(Parity::even);
              hipDeviceSynchronize();
            }
            //----------
            
             //--  onsites(ODD) {
            //--                  for (int i=0; i<N; i++)
            //--                          f_2[X] = f_2[X]*f_2[X + e_z]*f_2[X - e_z].dagger();
            //--              }
            {
              Field<SU<3, float>> & HILA_field_f_2 = f_2;
              HILA_field_f_2.check_alloc();
              if (HILA_field_f_2.fs->mylattice.ptr() != lattice.ptr()) {hila::out0 << "File /projappl/project_462000949/lumi-po-benchmark/benchmarks/matrix-test/src/matrix-test.cpp on line 39: Field f_2 initialized on different lattice\n";
                hila::terminate(1);
              }
              if (!HILA_field_f_2.is_initialized(ALL)) {
                hila::out0 << "File /projappl/project_462000949/lumi-po-benchmark/benchmarks/matrix-test/src/matrix-test.cpp on line 39:\n Variable f_2 is not properly initialized\n";
                hila::terminate(1);
              }
              dir_mask_t  _dir_mask_ = 0;
              roctxRangePush("Start pack");
              auto [var_1, mpi_1]= HILA_field_f_2.start_gather_test(e_z, Parity::odd);
              auto [var_2, mpi_2]= HILA_field_f_2.start_gather_test(- e_z, Parity::odd);
              gpuEventRecord(halo_event_2, halo_stream);
              roctxRangePop();
              roctxRangePush("Kernel bulk");
              const lattice_struct & loop_lattice = lattice.ref();
              int _hila_loop_begin = loop_lattice.loop_begin(Parity::odd);
              int _hila_loop_end = loop_lattice.loop_end(Parity::odd);
              int N_blocks = (_hila_loop_end - _hila_loop_begin + N_threads - 1)/N_threads;
              hipLaunchKernelGGL(HILA_kernel_main_K32, dim3(N_blocks), dim3(N_threads), 0, bulk_stream, _hila_loop_begin, _hila_loop_end, HILA_field_f_2.fs->payload);

              roctxRangePop();
              roctxRangePush("MPI send");
              gpuStreamWaitEvent(halo_stream, halo_event_2, 0);
              perform_mpi_send(mpi_1);
              perform_mpi_send(mpi_2);
              roctxRangePop();
              roctxRangePush("Wait and unpack");
              HILA_field_f_2.wait_gather(e_z, Parity::odd);
              HILA_field_f_2.wait_gather(- e_z, Parity::odd);
              roctxRangePop();
              check_device_error("HILA_kernel_main_K31");
              
              hila::set_allreduce(true);
              HILA_field_f_2.mark_changed(Parity::odd);
              hipDeviceSynchronize();
            }
            //----------
            
    }
    roctxRangePop();
    hila::synchronize();
    f_2_timer.stop();