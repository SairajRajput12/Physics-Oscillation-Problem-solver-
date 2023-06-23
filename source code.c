// here i have made this project using c language and using file handling concept of C and math library.


#include <stdio.h>
#include<math.h>
int
main ()
{
  int wish, choice;
  FILE *ptr;
  double ut, amplitude, t, t1, w, phase, velocity, dsr, u1;
  double ams, dms, acc, wms, phase2, phase3, vt;
  //Tushar variables  
  char st[33];
  int ch, def;
  double a, b, c, nt, M, v2, p, f;
  double pan, m, k, d, wn, T, h, v1, sq, n, mt;
  double ta, na, bha, ph;
  double Time_period, angular_velocity, uo, uo1;
  double tp;
  double zeta, c3, damp, c2, c4, c5, sum;
  double u0;
  double wd, u2, j, delta;
  double pr, tf, pf, tr, ui;

//sairaj 
  double ss, ti;
  double uj;
  //harsh
  double F_0, omega, A, omega_0, omega_r;
  ptr = fopen ("oscillation.txt", "w");
  printf
    ("welcome ti vibration guide here you will get to know about various types of vibration");
  printf
    ("for getting relationship between velocity,displacement and acceleration with time\n");
  printf ("press 0\n");
  printf ("  \n");
  printf ("for knowing about undamped oscillation enter 1\n");
  printf ("    \n");
  printf ("for getting knowladge about its application enter 2\n");
  printf ("     \n");
  printf ("for getting knowladge about damped oscillation enter 3\n");
  printf ("       \n");
  printf
    ("for getting knowladge about forced harmonic oscillation enter 4\n");

  scanf ("%d", &wish);

  switch (wish)
    {
    case 0:

      printf
	("entre the amplitude by which you are going know about vinlbration\n");
      scanf ("%lf %s", &amplitude, st);
      fprintf (ptr,
	       "The amplitude of a SHM can be defined as the maximum displacement of a particle from its mean position.\n");
      fprintf (ptr,
	       "you have entered the value of amplitude which is %f meter\n",
	       amplitude);
      /*writting defination with inputs */

      printf ("   \n");
      fprintf (ptr,
	       "The minimum time after which the particle keeps on repeating its motion is known as the time period\n");
      printf
	("entre the time period by which you are going know about vinlbration\n");
      scanf ("%lf %s", &t, st);
      printf ("   \n");
      fprintf (ptr,
	       "you have entered the value of time period which is %f seconds\n",
	       t);
      fprintf (ptr, ".        \n");
      w = 6.28 / t;
      fprintf (ptr,
	       "the rate of change of angular position of a rotating body about its center of rotation\n");
      fprintf (ptr, "    \n");
      printf ("thus the value of angular velocity is %f\n", w);
      fprintf (ptr,
	       "thus you have got the value of angular velocity which is:%f rad per second\n",
	       w);
      printf ("   \n");
      printf ("just entre what was initial displacement at t=0seconds");
      scanf ("%lf %s", &ut, st);
      fprintf (ptr,
	       "you have observed initial displacement of particle which is %f meter at time t=0 second\n",
	       ut);
      printf ("   \n");
      c3 = (sin (w * t1 - phase));
      ut = amplitude * c3;
      printf
	("enter time upto which you want to get displacement for finding phase of particle\n");
      scanf ("%lf %s", &t1, st);
      fprintf (ptr, "        \n");
      fprintf (ptr,
	       "phase tells us the condition of particle when observer started observing particle\n");
      fprintf (ptr, "        \n");
      phase = w * t1 - asin (ut / amplitude);
      printf ("thus the value of phase at any time t1=%f seconds  is  %f\n",
	      t1, phase);
      fprintf (ptr,
	       "thus the value of phase at any time t1=%f seconds  is  %f\n",
	       t1, phase);
      fprintf (ptr, "    \n");
      for (float i = 0; i <= t1; i++)
	{
	  ut = amplitude * sin (w * i - phase);
	  printf ("thus the value of displacement at time t=%f is is %f\n", i,
		  ut);
	  fprintf (ptr,
		   "thus the value of displacement at time t=%f is is %f\n",
		   i, ut);
	}
      fprintf (ptr, "    \n");
      printf ("ut=%f*sin(%f*t1-%f)\n", amplitude, w, phase);
      fprintf (ptr, "ut=%f*sin(%f*t1-%f)\n", amplitude, w, phase);
      //to set a programme for velocity.

      ams = pow (amplitude, 2);
      dms = pow (ut, 2);
      dsr = ams - dms;
      velocity = sqrt (dsr);
      fprintf (ptr, "    \n");
      fprintf (ptr,
	       "genrally the velocity is the rate of change in velocity per unit time\n");
      fprintf (ptr, " which can be given by dut/dt and we get\n");
      printf ("the equation for velocity is amplitude*w*cos(wt-phase)\n");
      fprintf (ptr, "    \n");
      fprintf (ptr, "in case of velocity amplitude*w will be max.velocity");
      fprintf (ptr, " which is %f\n", w * amplitude);
      fprintf (ptr, ".   \n");
      fprintf (ptr,
	       "the equation for velocity is amplitude*w*cos(wt-phase)\n");
      vt = amplitude * w * cos (w * t1 - phase);
      for (float i = 0; i <= t1; i++)
	{

	  vt = amplitude * w * cos (w * i - phase);
	  printf ("thus the value of velocity at time t=%f is is %f\n", i,
		  vt);

	}
      printf ("ut=%f*%f*cos(%f*t1-%f)\n", amplitude, w, w, phase);
      fprintf (ptr,
	       "thus i have concluded that the phase difference between velocity and displacement is 90degree\n");
      fprintf (ptr, "velocity leads displacement by 90 degree\n");
      phase2 = phase + 1.57;
      fprintf (ptr, "    \n");
      fprintf (ptr, "the phase of velocity of particle is %f\n ", phase2);
//to program for acceleration.//
      fprintf (ptr,
	       "the acceleration is rate of chamge of velocity w.r.t time \n");
      fprintf (ptr, "which can be given by dvt/dt\n");
      wms = pow (w, 2);
      acc = -1 * wms * ut;
      printf
	("the equation of displacement is -amplitude*(w) CB2*sin(wt-phase)\n");
      fprintf (ptr,
	       "the equation of acceleration is -amplitude*(w) CB2*sin(wt-phase)\n");
//printf("the value of acceleration at t=%f is %f",i,acc);  
      for (float i = 0; i <= t1; i++)
	{
	  acc = -1 * amplitude * (wms) * sin (w * i - phase);
	  printf ("thus the value of acceleration at time t=%f is is %f\n", i,
		  acc);
	  fprintf (ptr,
		   "thus the value of acceleration at time t=%f is is %f\n",
		   i, acc);

	}
      printf ("ut=-%f*%f*sin(%f*t1-%f)\n", amplitude, wms, w, phase);
      fprintf (ptr, "ut=-%f*%f*sin(%f*t1-%f)\n", amplitude, wms, w, phase);

      fprintf (ptr,
	       "thus i have concluded that the phase difference between acceleration and displacement is 180degree\n");
      fprintf (ptr, "acceleration leads displacement by 180degree\n");
      phase3 = phase + 3.14;
      fprintf (ptr, "the phase of velocity of acceleration is %f\n ", phase3);
      fprintf (ptr, "and also,acceleration leads velocity by 90 degree\n");

      fprintf (ptr,
	       "thus i have successfully learned the displacement, velocity and acceleration relation");
      printf ("by giving values from user, thank you");
      break;



    case 1:			//ammar khan
      fprintf (ptr,
	       "in undamped oscillation the resistive force are assumed to be zero \n");
      fprintf (ptr, "the amplitude in motion remains constant\n");
      printf
	("plz choose your desired way for angular frequency if you have calculated\n");
      printf
	("then press 1 otherwise you can interprete relation by giving time period\n");
      printf ("to the compiler\n");
      printf
	("if you entre 1 then it will ask you to enter angular velocity\n");
      printf ("otherwise it will ask you to enter time period\n");
      scanf ("%lf", &tp);
      if (tp == 1)
	{

	  printf ("entre the value of angular velocity\n");
	  scanf ("%lf %s", &angular_velocity, st);
	  fprintf (ptr,
		   "the rate of change of angular position of a rotating body about its center of rotation\n");
	  fprintf (ptr, "    \n");
	  Time_period = 6.28 / angular_velocity;
	  fprintf (ptr,
		   "The minimum time after which the particle keeps on repeating its motion is known as the time period\n");

	  printf ("the time period of particle is %lf sec\n", Time_period);
	  fprintf (ptr,
		   "you have entred value of angular velocity is %f radper second\n",
		   angular_velocity);
	  fprintf (ptr, "the time period of particle is %f sec\n",
		   Time_period);

	}
      else
	{
	  fprintf (ptr,
		   "The minimum time after which the particle keeps on repeating its motion is known as the time period\n");

	  printf ("entre time period\n");
	  scanf ("%lf %s", &Time_period, st);
	  f = 1 / Time_period;
	  angular_velocity = 6.28 * f;
	  fprintf (ptr,
		   "the rate of change of angular position of a rotating body about its center of rotation\n");
	  fprintf (ptr, "    \n");
	  printf ("the value of angular velocity is-%f rad per sec\n",
		  angular_velocity);
	  fprintf (ptr, "you have entred value of  time period is %f\n",
		   Time_period);
	  fprintf (ptr, "the value of angular velocity is:%f rad per sec\n",
		   angular_velocity);
	}
      t1 = uo1 / angular_velocity;
      printf ("entre value of displacement at t=0sec\n");
      scanf ("%lf %s", &uo, st);
      printf ("   \n");
      printf ("entre value of velocity at t=0 sec\n");
      scanf ("%lf %s", &uo1, st);
      printf ("    \n");


      t1 = uo1 / angular_velocity;
      phase = atan (uo / t1);
      c = pow (t1, 2);
      b = pow (uo, 2);
      a = b + c;
      amplitude = sqrt (a);
      fprintf (ptr,
	       "The amplitude of a SHM can be defined as the maximum displacement of a particle from its mean position.\n");

      printf ("you have entred value of amplitude is %f metre \n", amplitude);
      fprintf (ptr, "you have entred value of amplitude is %f metre \n",
	       amplitude);

      printf
	("entre value of time upto which you want to get displacement\n");
      scanf ("%lf %s", &t, st);

      d =
	uo * sin (angular_velocity * t) +
	(uo1 / angular_velocity) * cos (angular_velocity * t);
      for (int i = 0; i <= t; i++)
	{

	  d = amplitude * (sin (angular_velocity * i - phase));

	  fprintf (ptr, "the displacement of particle is %f meters\n", d);

	}
      fprintf (ptr,
	       "thus from reading i have concluded that the amplitude of particle remains constant");
      fprintf (ptr, "due to absence of Resistive forces");
      fprintf (ptr,
	       "thus we have got the equation of particle is %f*sin(%f*t-%f)\n",
	       amplitude, angular_velocity, phase);

      break;

    case 2:			//atharva walkhede
      printf
	("you have entered '0'which means you have selected to present vibration example in spring mass problem \n");
      printf ("      \n");
      printf
	("you have entered '1'which means you have selected to present vibration example in bullet and block problem \n");
      printf ("      \n");

      printf ("enter your choice\n");
      scanf ("%d", &def);
      switch (def)
	{
	case 0:
	  printf
	    ("you have entered '0'which means you have selected to present vibration example in spring mass problem \n");
	  printf ("plz enter value of mass\n");
	  scanf ("%lf %s", &m, st);
	  printf ("      \n");
	  printf ("plz enter value of spring constant\n");
	  scanf ("%lf %s", &k, st);

	  d = k / m;
	  wn = sqrt (d);
	  fprintf (ptr,
		   "the angular velocity is rate of change in angular position\n");
	  fprintf (ptr, ".   \n");
	  printf
	    ("thus,the value of angular velocity is: %f radian per second\n",
	     wn);
	  printf ("      \n");
	  fprintf (ptr,
		   "thus,the value of angular velocity is: %f radian per second\n",
		   wn);
	  fprintf (ptr, ".   \n");
	  T = wn / 6.28;
	  printf ("thus,the value of time period is: %f  second\n", T);
	  fprintf (ptr, "thus,the value of time period is: %f  second\n", T);

	  printf
	    ("if you entre 1 then you will find the velocity by using momentum conservation\n");
	  printf ("      \n");
	  printf
	    ("if you entre 2 then you will find the velocity by using energy conservation\n");

	  printf ("enter your choice\n");
	  scanf ("%d", &ch);

	  if (nt == 1)
	    {
	      printf ("entre the mass of bigger object\n");
	      scanf ("%lf %s", &m, st);
	      printf ("      \n");
	      printf ("entre the mass of smaller object\n");
	      scanf ("%lf %s", &M, st);
	      printf ("      \n");
	      mt = m + M;
	      v1 = (m * v2) / mt;
	      printf ("entre the initial velocity of smaller object\n");
	      scanf ("%lf %s", &v2, st);
	      printf ("      \n");
	      fprintf (ptr,
		       "the final velocity of other body is give. by  v1=(m*v2)/mt");

	      printf ("the value of  velocity is: %lf  meter per second\n",
		      v1);


	    }

	  else
	    {
	      printf ("      \n");
	      printf ("plz enter value of distance traveled by particle\n");
	      scanf ("%lf %s", &h, st);
	      n = 2 * (9.8) * h;
	      v1 = sqrt (n);
	      printf ("velocity of particle is %f m/sec\n", v1);
	      printf ("      \n");


	      printf ("displacement of particle at t=0 sec will be\n");
	      printf ("      \n");
	      printf ("u(t)=u(o)*sin(wt)+(v1/wn)cos(wt) is equation\n");
	      printf ("      \n");
	      printf ("at t=0 ,u(t) will be zero, therefore\n");
	      printf ("      \n");
	      t = 1;
	      d = k / m;
	      wn = sqrt (d);
	      bha = pow (cos (wn * t), 2);
	      printf ("bha %f", bha);
	      na = 1 / bha;
	      printf ("na -%f", na);
	      ta = na - 1;

	      pan = sqrt (ta);
	      uo = -1 * (v1 / wn) * pan;
	      printf ("displacement of particle at t=0sec is %f m/sec\n", uo);
	      t1 = v1 / wn;
	      c = pow (t1, 2);
	      b = pow (uo, 2);
	      a = b + c;
	      amplitude = sqrt (a);
	      printf ("amplitude of particle is %lf m\n", amplitude);
	      printf ("      \n");
	      phase = atan (uo / t1);
	      printf ("phase of particle is %f m\n", phase);
	      printf ("      \n");
	      printf ("enter your time\n");
	      scanf ("%lf %s", &t, st);
	      printf ("      \n");
	    }

	  for (float i = 0; i <= t; i++)
	    {
	      ut = amplitude * sin (wn * i - phase);

	      printf ("the displacement at t=%f is %f\n", i, ut);

	    }
	  printf ("the equation of motion is %f*sin(%f*t-%f)\n", amplitude,
		  wn, phase);
	  break;


	case 1:
	  printf ("entre the mass of bullete object\n");
	  scanf ("%lf %s", &m, st);
	  printf ("      \n");
	  printf ("entre the mass of block object\n");
	  scanf ("%lf %s", &M, st);
	  printf ("      \n");
	  printf ("entre the velocity of bullet fired object\n");
	  scanf ("%lf %s", &v2, st);
	  printf ("      \n");
	  mt = m + M;
	  v1 = (m * v2) / mt;
	  fprintf (ptr, "the combine velocity will be v1=(m*v2)/mt\n");

	  printf ("plz enter value of spring constant\n");
	  scanf ("%lf %s", &k, st);
	  printf ("      \n");
	  printf
	    ("we know that,the initial displacement and velocity of block is zero\n");
	  printf ("because it is on rest\n");
	  printf ("therefore the phase of particle initial is zero\n");
	  printf ("      \n");
	  printf ("the value of  velocity is: %f  second\n", v1);
	  printf ("      \n");
//  printf("the equation of motion is %f*sin(%f*t)\n", amplitude,wn);
	  p = k / mt;
	  wn = sqrt (p);
	  printf ("      \n");
	  printf ("the value of angular velocity is: %f  second\n", wn);
	  f = 6.28 / wn;
	  T = 1 / f;
	  printf ("      \n");
	  printf
	    ("the value of frequency  is: %f persecond and time period is %f\n",
	     wn, T);
	  printf ("      \n");
	  amplitude = v1 / wn;

	  for (float i = 0; i <= t; i++)
	    {
	      ut = amplitude * sin (wn * t);

	      printf ("the displacement at t=%f is %lf\n", i, ut);

	    }
	  printf ("      \n");
	  printf ("the equation of motion is %f*sin(%f*t)\n", amplitude, wn);
	  break;

	default:
	  printf ("bye\n");

	}
      break;

    case 3:			//sairaj Rajput
      {
	printf ("for finding value of damped angular velocity entre 0\n");

	printf ("     \n");
	printf ("for finding value of amplitude entre 1\n");

	printf ("     \n");

	printf ("for finding value of phase difference entre 2\n");
	//printf ("for getting angular velocity enter 3\n"); 
	printf ("     \n");

	printf ("for getting value of damping ratio just entre 3\n");

	printf ("   \n");
	printf ("for getting damping dactor just entre 4\n");

	printf ("      \n");
	printf ("entre your choice of topic\n");


	printf ("      \n");

	scanf ("%d", &choice);

	printf ("     \n");


	printf ("    \n");
	printf ("    \n");

	switch (choice)
	  {

	  case 0:


	    printf
	      ("you have selected to find value of damped angular velocity\n ");
	    printf ("    \n");
	    printf ("for finding damped angular velocity you have 2 ways\n");
	    printf ("     \n");

	    printf
	      ("on entering '1' you have to give  parameters like stiffness and mass of block for finding angular velocity \n");
	    printf ("      \n");
	    printf ("otherwise it will ask you value of frequency\n");
	    printf ("    \n");
	    scanf ("%lf", &d);

	    {
	      if (d == 1)
		{
		  printf
		    ("now you have to find value of angular velocity by experiment");
		  printf ("   \n");
		  printf ("entre value of mass of body");
		  printf ("    \n");
		  scanf ("%lf %s", &m, st);
		  printf ("entre value of stiffness constant of spring");
		  printf ("    \n");
		  scanf ("%lf %s", &k, st);

		  ti = k / m;
		  w = pow (ti, 0.5);
		  f = w / 6.28;
		  T = 1 / f;
		  printf ("    \n");
		  fprintf (ptr,
			   "you have got the values of time period:-%f seconds \n",
			   T);
		  printf ("    \n");
		  fprintf (ptr,
			   "you have got the values of frequency:-%f per second\n",
			   f);
		  printf ("     \n");
		  fprintf (ptr,
			   "you have got the values of angular velocity:-%f radian per second\n",
			   w);
		  printf ("    \n");
		}

	      else
		{
		  printf ("entre value of time period :\n");

		  scanf ("%lf %s", &T, st);
		  printf ("     \n");
		  printf ("the value of frequency is %f\n", 1 / T);
		  w = 6.28 * (1 / T);
		  printf ("    \n");
		  fprintf (ptr, "the value angular velocity is:%f\n", w);
		}
	    }
	    printf ("entee zeta\n");
	    scanf ("%lf %s", &zeta, st);
	    if (zeta == 1)
	      {
		fprintf (ptr, "there is no oscillation\n");

		fprintf (ptr, "    \n");
		fprintf (ptr, "time period is infinite\n");
	      }
	    else if (zeta > 1)
	      {
		fprintf (ptr, "it is overdamped case \n");
	      }
	    else
	      {
		fprintf (ptr, "it is underdamped case\n");

		c3 = pow (zeta, 2);
		d = (1 - c3);
		wd = (pow (d, 0.5));
		damp = w * wd;

		fprintf (ptr, "   \n");
		fprintf (ptr,
			 "you have got the damped angular frequency is \n:%f",
			 damp);

	      }

	    break;

	  case 1:

	    printf (" you have choiced to find the value of amplitude\n");
	    printf ("   \n");
	    printf ("for finding amplitude  you have 2 ways\n");
	    printf ("   \n");
	    printf
	      ("on entering '1' you have to give  parameters like stiffness and mass of block for finding angular velocity \n");
	    printf ("    \n");
	    printf ("otherwise it will ask you value of frequency\n");
	    scanf ("%lf %s", &d, st);
	    {

	      if (d == 1)
		{
		  printf
		    ("now you have to find value of angular velocity by experiment\n");
		  printf ("    \n");
		  printf ("entre value of mass of body\n");
		  scanf ("%lf", &m);
		  printf ("    \n");
		  printf ("entre value of stiffness constant of spring");
		  scanf ("%lf %s", &k, st);
		  ti = k / m;
		  w = pow (ti, 0.5);
		  f = w / 6.28;
		  T = 1 / f;

		  printf ("time period you got-%fseconds", T);
		  fprintf (ptr,
			   "you have got the values of time period:-%f seconds \n",
			   T);
		  fprintf (ptr, "   \n");
		  fprintf (ptr,
			   "you have got the values of frequency:-%f per seconds\n",
			   f);
		  fprintf (ptr, "   \n");
		  printf ("   \n");
		  fprintf (ptr,
			   "you have got the values of angular velocity:-%f radper seconds\n",
			   w);
		  fprintf (ptr, "   \n");
		}


	      else
		{
		  printf ("entre value of time period:\n");
		  scanf ("%lf %s", &T, st);
		  f = 1 / T;
		  w = 6.28 * f;

		  fprintf (ptr,
			   "the value angular velocity is:%f rad per second\n",
			   w);
		}


	      printf ("entre diaplacement at t=0\n");
	      scanf ("%lf %s", &u1, st);
	      printf ("   \n");
	      printf ("entee velocity at t=0 \n");
	      scanf ("%lf %s", &u2, st);
	      printf ("    \n");
	      printf ("entee zeta\n");
	      scanf ("%lf", &zeta);
	      printf ("    \n");
	      if (zeta == 1)
		{

		  fprintf (ptr, "there is no damped oscillation\n");
		  fprintf (ptr, "time period is infinite\n");

		}
	      else if (zeta > 1)
		{
		  fprintf (ptr, "it is overdamped case \n");
		}
	      else
		{
		  fprintf (ptr, "it is underdamped case\n");

		  c3 = pow (zeta, 2);
		  d = (1 - c3);
		  wd = (pow (d, 0.5));
		  damp = w * wd;
		  c4 = pow (u1, 2);
		  n = ((u1 + zeta * u2 * w) / damp);

		  c5 = pow (n, 2);
		  sum = c5 + c4;
		  /*"I HAVE CONCLUDED THAT WHEN ZETA >0 MEANS WE GET IMAGINARY NUMBER
		     COMPILER IS UNABLE TO DEFINE IT"/* " */
		  u0 = pow (sum, 0.5);
		  fprintf (ptr, "the value of amplitude is :%f\n", u0);
		  printf ("    \n");
		  printf (" entre time\n");
		  scanf ("%lf %s", &t, st);

		  for (float i = 0.5; i <= t; i++)
		    {
		      ph = 0;
		      t1 = exp (zeta * i * w);
		      ui = u0 * sin (w * i - ph) * t1;
		      fprintf (ptr,
			       "the displacement at time t=%f seconds is %f m\n",
			       i, ui);
		    }

		}
	      break;

	  case 2:
	      printf
		("you have choiced  to find the value of phase difference of\n ");
	      printf ("    \n");
	      printf ("for finding phase difference   you have 2 ways\n");
	      printf ("    \n");
	      printf
		("on entering '1' you have to give  parameters like stiffness and mass of block for finding angular velocity \n");
	      printf ("   \n");
	      printf ("otherwise it will ask you value of frequency\n");
	      scanf ("%lf", &d);
	      printf ("   \n");

	      {

		if (d == 1)
		  {
		    printf
		      ("now you have to find value of angular velocity by experiment");
		    printf ("     \n");
		    printf ("entre value of mass of body");
		    scanf ("%lf %s", &m, st);
		    printf ("    \n");
		    printf ("entre value of stiffness constant of spring");
		    scanf ("%lf %s", &k, st);
		    ti = k / m;
		    w = pow (ti, 0.5);
		    f = w / 6.28;
		    T = 1 / f;
		    fprintf (ptr,
			     "you have got the values of time period:-%f\n",
			     T);
		    printf ("     \n");
		    fprintf (ptr,
			     "you have got the values of frequency:-%f\n", f);
		    printf ("    \n");
		    fprintf (ptr,
			     "you have got the values of angular velocity:-%f\n",
			     w);
		    printf ("    \n");

		  }

		printf ("entre diaplacement at t=0\n");
		scanf ("%lf %s", &u1, st);
		printf ("   \n");
		printf ("entee velocity at t=0 \n");
		scanf ("%lf %s", &u2, st);
		printf ("    \n");
		printf ("entee zeta\n");
		scanf ("%lf %s", &zeta, st);
		printf ("    \n ");
		if (zeta == 1)
		  {

		    fprintf (ptr, "there is no oscillation\n");
		    fprintf (ptr, "   \n");
		    fprintf (ptr, "time period is infinite\n");

		  }
		else if (zeta > 1)
		  {
		    fprintf (ptr, "it is overdamped case \n");
		  }
		else
		  {
		    fprintf (ptr, "it is underdamped case\n");
		    printf ("entre value of frequency");
		    scanf ("%lf %s", &f, st);
		    w = 6.28 * f;
		    c3 = pow (zeta, 2);
		    d = (1 - c3);
		    wd = (pow (d, 0.5));
		    damp = w * wd;
		    fprintf (ptr,
			     "the value of damped angular velocity is \n:%f",
			     damp);
		    c4 = pow (u1, 2);
		    n = ((u1 + zeta * u2 * w) / wd);
		    ti = -1 * (u1 / n);
		    ss = atan (ti);
		    fprintf (ptr, "the phase difference is %f\n", ss);


		  }
	      }

	      break;
	  case 3:

	      printf ("you have choiced to find value of damping ratio\n");
	      printf ("    \n");
	      printf ("entre the value of jth amplitude\n ");
	      scanf ("%lf %s", &uj, st);
	      printf ("    \n");
	      printf ("entre value of j\n");
	      scanf ("%lf %s", &j, st);
	      printf ("    \n");
	      printf ("entre value of ith amplitude \n");
	      scanf ("%lf %s", &u1, st);

	      fprintf (ptr, "the damping ratio is %f\n", delta);
	      ph = u1 / uj;
	      delta = (1 / j) * 2.303 * log (ph);
	      c3 = pow (delta, 2);
	      c2 = pow (39.43, 2);
	      tr = c2 + c3;
	      pf = pow (tr, 0.5);
	      zeta = delta / pf;

	      fprintf (ptr,
		       "you have got the value of damping factor which is %f\n",
		       zeta);

	      break;

	  default:
	      {


	      }
	    }
	  }

	break;

    case 4:			//harshwardhan kulkarni
	printf ("Forced Harmonic Oscillation\n");

	printf ("Enter mass:");
	scanf ("%lf %s", &m, st);

	printf ("Enter force constant:");
	scanf ("%lf %s", &k, st);

	printf ("Enter damping coefficient:");
	scanf ("%lf %s", &b, st);

	printf ("Enter angular frequency of the external force:");
	scanf ("%lf %s", &omega, st);

	printf ("Enter Force per unit mass:");
	scanf ("%lf %s", &F_0, st);

	omega_0 = sqrt (k / m);
	//Natural Angular frequency

	A =
	  (F_0 / m) /
	  (sqrt
	   ((omega * omega - omega_0 * omega_0) * (omega * omega -
						   omega_0 * omega_0) +
	    ((b * omega) / m) * ((b * omega) / m)));
	//Amplitude

	omega_r = sqrt (omega_0 * omega_0 - ((b * b) / (4 * m * m)));
	//Resonant frequency

	printf ("The amplitude of the oscillation is %f\n", A);
	printf ("The natural angular frequency of the oscillation is %f\n",
		omega_0);
	printf ("The resonant frequency of the oscillation is %f\n\n",
		omega_r);
	printf
	  ("*These calculations are valid for an external force of the form F=F_0.sin(N)t)");


	break;
    default:
    //printf("congrates you got knowladge about oscillation"); 
      }
    }
fclose(ptr);
return 0; 
}
