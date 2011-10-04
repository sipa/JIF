    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : q*2.5;
    double qq = argc > 4 ? strtod(argv[4],NULL) : q*3;
    double factor = argc > 5 ? strtod(argv[5],NULL) : 1.0*(2*q+qi+qq)/8.0;
    double factor_s = argc > 6 ? strtod(argv[6],NULL) : factor/5.0;
    double qzy = clamp(argc > 7 ? strtod(argv[7],NULL) : 1+(q+1)/3,1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 8 ? strtod(argv[8],NULL) : 1+(qi+1)/3,1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 9 ? strtod(argv[9],NULL) : 1+(qq+1)/3,1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 10 ? strtol(argv[10],NULL,10) : clamp(9 - qzy/2,1,8),8,8);  // using all bits always, otherwise bugs in symbol out
    sb[1] = clamp(argc > 11? strtol(argv[11],NULL,10) : clamp(10 - qzi/2,1,9),9,9);
    sb[2] = clamp(argc > 12? strtol(argv[12],NULL,10) : clamp(10 - qzq/2,1,9),9,9);

/*
    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : q*2.5;
    double qq = argc > 4 ? strtod(argv[4],NULL) : q*3;
    double factor = argc > 5 ? strtod(argv[5],NULL) : 1.0*(2*q+qi+qq)/6.0;
    double factor_s = argc > 6 ? strtod(argv[6],NULL) : factor/5.0;
    double qzy = clamp(argc > 7 ? strtod(argv[7],NULL) : 1+(q+1)/3,1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 8 ? strtod(argv[8],NULL) : 1+(qi+1)/3,1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 9 ? strtod(argv[9],NULL) : 1+(qq+1)/3,1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 10 ? strtol(argv[10],NULL,10) : clamp(9 - qzy/2,1,8),0,8);
    sb[1] = clamp(argc > 11? strtol(argv[11],NULL,10) : clamp(10 - qzi/2,1,9),0,9);
    sb[2] = clamp(argc > 12? strtol(argv[12],NULL,10) : clamp(10 - qzq/2,1,9),0,9);

    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : q*2.5;
    double qq = argc > 4 ? strtod(argv[4],NULL) : q*3;
    double factor = argc > 5 ? strtod(argv[5],NULL) : 1.0*(2*q+qi+qq)/4.0;
    double factor_s = argc > 6 ? strtod(argv[6],NULL) : factor/5.0;
    double qzy = clamp(argc > 7 ? strtod(argv[7],NULL) : 1+(q+1)/3,1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 8 ? strtod(argv[8],NULL) : 1+(qi+1)/3,1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 9 ? strtod(argv[9],NULL) : 1+(qq+1)/3,1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 10 ? strtol(argv[10],NULL,10) : clamp(9 - qzy/2,1,8),0,8);
    sb[1] = clamp(argc > 11? strtol(argv[11],NULL,10) : clamp(10 - qzi/2,1,9),0,9);
    sb[2] = clamp(argc > 12? strtol(argv[12],NULL,10) : clamp(10 - qzq/2,1,9),0,9);


    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : q*2.5;
    double qq = argc > 4 ? strtod(argv[4],NULL) : q*3;
    double factor = argc > 5 ? strtod(argv[5],NULL) : 1.0*(2*q+qi+qq)/3.0;
    double factor_s = argc > 6 ? strtod(argv[6],NULL) : factor/5.0;
    double qzy = clamp(argc > 7 ? strtod(argv[7],NULL) : (q+2)/3,1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 8 ? strtod(argv[8],NULL) : (qi+2)/3,1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 9 ? strtod(argv[9],NULL) : (qq+2)/3,1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 10 ? strtol(argv[10],NULL,10) : clamp(9 - qzy,2,8),0,8);
    sb[1] = clamp(argc > 11? strtol(argv[11],NULL,10) : clamp(10 - qzi,2,9),0,9);
    sb[2] = clamp(argc > 12? strtol(argv[12],NULL,10) : clamp(10 - qzq,2,9),0,9);


    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : q*2.5;
    double qq = argc > 4 ? strtod(argv[4],NULL) : q*3;
    double factor = argc > 5 ? strtod(argv[5],NULL) : 2.0*(2*q+qi+qq)/3.0;
    double factor_s = argc > 6 ? strtod(argv[6],NULL) : factor/10.0;
    double qzy = clamp(argc > 7 ? strtod(argv[7],NULL) : (q+2)/4,1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 8 ? strtod(argv[8],NULL) : (qi+2)/4,1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 9 ? strtod(argv[9],NULL) : (qq+2)/4,1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 10 ? strtol(argv[10],NULL,10) : clamp(9 - qzy,2,8),0,8);
    sb[1] = clamp(argc > 11? strtol(argv[11],NULL,10) : clamp(10 - qzi,2,9),0,9);
    sb[2] = clamp(argc > 12? strtol(argv[12],NULL,10) : clamp(10 - qzq,2,9),0,9);
*/

/*
        // orig:
    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : q*3;
    double qq = argc > 4 ? strtod(argv[4],NULL) : q*4;
    double factor = argc > 5 ? strtod(argv[5],NULL) : (q+qi+qq)/2;
    double qzy = clamp(argc > 6 ? strtod(argv[6],NULL) : sqrt(q),1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 7 ? strtod(argv[7],NULL) : sqrt(qi),1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 8 ? strtod(argv[8],NULL) : sqrt(qq),1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 9 ? strtol(argv[9],NULL,10) : (9 - qzy),0,8);
    sb[1] = clamp(argc > 10? strtol(argv[10],NULL,10) : (10 - qzi),0,9);
    sb[2] = clamp(argc > 11? strtol(argv[11],NULL,10) : (10 - qzq),0,9);

        // 2:
    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : q*3;
    double qq = argc > 4 ? strtod(argv[4],NULL) : q*4;
    double factor = argc > 5 ? strtod(argv[5],NULL) : sqrt(q+qi+qq);
    double qzy = clamp(argc > 6 ? strtod(argv[6],NULL) : q/3,1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 7 ? strtod(argv[7],NULL) : qi/3,1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 8 ? strtod(argv[8],NULL) : qq/3,1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 9 ? strtol(argv[9],NULL,10) : clamp(9 - qzy,2,8),0,8);
    sb[1] = clamp(argc > 10? strtol(argv[10],NULL,10) : clamp(10 - qzi,2,9),0,9);
    sb[2] = clamp(argc > 11? strtol(argv[11],NULL,10) : clamp(10 - qzq,2,9),0,9);

        // 4:
    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : q*3;
    double qq = argc > 4 ? strtod(argv[4],NULL) : q*4;
    double factor = argc > 5 ? strtod(argv[5],NULL) : 1+(2*q+qi+qq)/2;
    double qzy = clamp(argc > 6 ? strtod(argv[6],NULL) : (q+2)/4,1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 7 ? strtod(argv[7],NULL) : (qi+2)/4,1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 8 ? strtod(argv[8],NULL) : (qq+2)/4,1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 9 ? strtol(argv[9],NULL,10) : clamp(9 - qzy,2,8),0,8);
    sb[1] = clamp(argc > 10? strtol(argv[10],NULL,10) : clamp(10 - qzi,2,9),0,9);
    sb[2] = clamp(argc > 11? strtol(argv[11],NULL,10) : clamp(10 - qzq,2,9),0,9);

*/
