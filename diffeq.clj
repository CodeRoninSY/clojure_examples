;; -----------------------------------------------------------------------------
;; diffeq.core
;;
;; ODE numerical methods
;;     :: Personal project for learning Clojure numerical modeling capabilities & shortbacks
;;
;; Author = SY
;; Date = <2017-04-29>
;; Last update = <2017-05-06>
;;
;; Methods Implemented:
;;   1. Euler
;;   2. Heun
;;   3. Adams-Bashfort 2-step
;;   4. Adams-Bashfort-Moulton 4th order predictor-corrector
;;   5. Runge-Kutta 2nd order A
;;   6. Runge-Kutta 2nd order B
;;   7. Runge-Kutta 4th order
;;   8. Runge-Kutta 4th order with error estimate
;;   9. Runge-Kutta-Fehlberg

;; -----------------------------------------------------------------------------

;; namespace for DiffEq
(ns diffeq.core
  (:gen-class)
  (:require [clojure.java.io :as io]
            [clojure.string :as strn]))

(def tofile "diffeq_out.dat")

;; Write the solution numerics to file
;;------------------------------------------------------------------------------
;; function "write-to" : write all generated f(x,t) values to file
(defn write-to [fname f header]
  "write to file function
    inputs:
      fname = file name to be written to
      f = function results to write
      header = header to write"
  (with-open [wrt (io/writer fname :append true)]
    (binding [*out* wrt]
      (println header)
      (doseq [q f]
        (println (apply format "t= %16.12f \t f= %16.12f \t ex= %16.12f \t err= %-8.3e" q))))))

;;------------------------------------------------------------------------------

;; Functions and exact solutions to funcs - for numerics comparison
;; x' = (3t - 2x)/t; x(1)=0
;;(defn f [x t]
;;(/ (- (* 3 t) (* 2 x)) t))

;; exact solution to x
;;(defn exact [t]
;;(/ (- (* t t t) 1) (* t t)))

;; x'= x sin(t); x(0)=-1.0
(defn ^double f [^double x ^double t]
  (* x (Math/sin t)))

(defn ^double exact [^double t]
  (* -1.0 (Math/exp (- 1 (Math/cos t)))))


;;------------------------------------------------------------------------------

;; parameter settings
(def a 0.0)  ;; t_start
(def b 10.0)  ;; t_end
(def t0 0.0) ;; t0
(def x0 -1.0) ;; x0

(def n 101)  ;; total no of steps

(def step (/ (- b a) n))   ;; step size

;;------------------------------------------------------------------------------
;; Euler method
(defn euler [f fx x0 a b h]
  "Euler's Method.
  Approximates y(time) in y'(time)=f(time,y) with y(a)=y0 and t=a..b and the step size h.
  Usage:
        (euler f fexact x0 a b h)
  Input:
    f : function of x and t f(x,t) equal to dx/dt
    x0 : x_start
    a : t0
    b : t_end
    h : step size
  "
  (println "\n\n*** EULER ***\n\n")
  (loop [t a
         x x0
         result []]
    (if (<= t b)
      (recur (+ t h) (+ x (* (f (+ t h) x) h))
             (conj result [(double t) (double x) (double (fx t)) (- x (fx t))]))
      result)))

;;------------------------------------------------------------------------------


;; Heun method
(defn heun [f fx x0 t0 step n]
  " Heun's method to solve y'=f(y,x) with y(x[0])=y0
    Usage: (heun x0 t0 step n)
    Input:
      f : function of x and t equal to dx/dt
      x0 : x_start
      t0 : t_start
      step: step value
      n : no of steps"
  (println "\n\n*** HEUN ***\n\n")
  (loop [x x0
         t t0
         i n
         result []]
    (println "step= " (- n i) ", t= " t ", x= " x ", exact= " (fx t))
    (if (>= i 0)
      (let [h step
            k1 (f x t)
            t1 (+ t (* (/ 2 3) h))
            x1 (+ x (* (/ 2 3) h k1))
            k2 (f x1 t1)
            ex (fx t)
            err (- x ex)]
        (recur (+ x (* 0.25 h (+ k1 (* 3 k2)))) (+ t h) (dec i)
               (conj result [t x ex err])))
      result)))




;;------------------------------------------------------------------------------

;; 2-step Adams-Bashford method
(defn adbash-2 [f fx x0 t0 step n]
  " two step adams-bashford method"

  (println "\n\n*** ADBASH-2 ***\n\n")
  (loop [
         j 0
         x x0
         t t0
         result []
         ]
    (println "step= " j ", t= " t ", x= " x ", exact= " (fx t))
    (if (<= j n)
      (let [
            h step
            f0 (f x t)
            k1 (* h f0)
            k2 (* h (f (+ x (* 0.5 k1)) (+ t (* 0.5 h))))
            k3 (* h (f (+ x (* 0.5 k2)) (+ t (* 0.5 h))))
            k4 (* h (f (+ x k3) (+ t h)))
            ex (fx t)
            ]
        ;(println "-- step = " j " , t= " t " , x= " x)
        (recur
          (inc j)
          (+ x (/ (+ k1 (* 2.0 (+ k2 k3)) k4) 6.0))
          (+ t h)
          (conj result [t x ex (- x ex)])
          )
        )
      result)))



;; -----------------------------------------------------------------------------
;; Adams-Bashfort-Moulton 4th order predictor-corrector method
(defn abm4-pc [f fx x0 t0 step n]
  "Adams-Bashforth-Moulton 4th order predictor-corrector method

    USAGE:
        (abm4-pc f f_ex x0 t0 step n)

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a  array.  In this
                case f must return a  array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or  array
                if a system of equations is being solved.
        t     - list or  array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - variable containing solution values corresponding to each
                entry in t var.

    NOTES:
        This function used the Adams-Bashforth-Moulton predictor-corrector
        method to solve the initial value problem

            dx
            -- = f(x,t),     x(t(1)) = x0
            dt

        at the t values stored in the t array (so the interval of solution is
        [t[0], t[N-1]].  The 4th-order Runge-Kutta method is used to generate
        the first three values of the solution.  Notice that it works equally
        well for scalar functions f(x,t) (in the case of a single 1st order
        ODE) or for vector functions f(x,t) (in the case of multiple 1st order
        ODEs).

    "

  (println "\n\n*** ABM4-PC ***\n\n")
  (loop [
         j 0
         x x0
         t t0
         result []  ;; results array
         ff ()  ;; store (f1, f2, f3) for using in predictor-corrector block
         ]
    ;(println  "-> Loop root :: step = " j " , t= " t " , x= " x ";  f1= " f1 " , f2= " f2 " , f3= " f3 )
    (cond
      ;; this cond case block utilizes the Runge-Kutta 4th for first 3 steps to
      ;; evaluate initial values for x
      (<= j (min 3 (- n 1))) (let [
                                   h step
                                   f0 (f x t)
                                   k1 (* h f0)
                                   k2 (* h (f (+ x (* 0.5 k1)) (+ t (* 0.5 h))))
                                   k3 (* h (f (+ x (* 0.5 k2)) (+ t (* 0.5 h))))
                                   k4 (* h (f (+ x k3) (+ t h)))
                                   ex (fx t)
                                   f1 f0
                                   f2 (first ff)
                                   f3 (second ff)
                                   ff (take 3 (conj ff f3 f2 f1))
                                   ]
                               (println  "-> RungeKutta 1st step :: step = " j " , t= " t " , x= " x
                                         ";; f0= " f0 ", f1= " f1 " , f2= " f2 " , f3= " f3)
                               (println ff)
                               (recur
                                 (inc j)
                                 (+ x (/ (+ k1 (* 2.0 (+ k2 k3)) k4) 6.0))
                                 (+ t h)
                                 (conj result [t x ex (- x ex)])
                                 (take 3 (conj ff f3 f2 f1))

                                 )
                               )
      ;; this cond case block utilizes the predictor-corrector steps after first
      ;; runge-kutta values
      (<= j n) (let [
                     h step
                     f0 (f x t)
                     f1 (first ff)
                     f2 (second ff)
                     f3 (last ff)
                     w (+ x (/ (* h (+ (* 55.0 f0) (* -59.0 f1)
                                                   (* 37.0 f2) (* -9.0 f3))) 24.0))
                     fw (f w (+ t h))
                     x1 (/ (* h (+ (* 9.0 fw) (* 19.0 f0) (* -5.0 f1) f2) ) 24.0)
                     f1 f0
                     f2 (first ff)
                     f3 (second ff)
                     ff (take 3 (conj ff f3 f2 f1))
                     ex (fx t)
                     ]
                 (println  "-> Adams-Bashfort-Moulton 2nd :: step = " j " , t= " t " , x= " x
                           ";; f0= " f0 ";  f1= " f1 " , f2= " f2 " , f3= " f3)
                 (println ff)
                 (recur
                   (inc j)
                   (+ x x1)
                   (+ t h)
                   (conj result [t x ex (- x ex)])
                   (take 3 ff)

                   )
                 )
      :else result)
    )
  )

;;------------------------------------------------------------------------------
;; Runge-Kutta 2nd order method
(defn rk2a [f fx x0 t0 step n]
  " RK2A - second order runge-kutta method to solve x'=f(x,t), x(t[0])=x0"
  (println "\n\n*** RUNGE-KUTTA 2A ***\n\n")
  (loop [x x0
         t t0
         i n
         result []]
    (println "step= " (- n i) ", t= " t ", x= " x " , exact= " (fx t))
    (if (>= i 0)
      (let [h step
            k1 (* h (/ (f x t) 2))
            ex (fx t)]
        (recur (+ x (* h (f (+ x k1) (+ t (/ h 2))))) (+ t h) (dec i)
               (conj result [(double t) (double x) (double ex) (- x ex)])))
      result)))

;;------------------------------------------------------------------------------
;; Runge-Kutta 2nd order method
(defn rk2b [f fx x0 t0 step n]
  " RK2B - second order runge-kutta method"
  (println "\n\n*** RUNGE-KUTTA 2B ***\n\n")
  (loop [x x0
         t t0
         i n
         result []]
    (println "step= " (- n i) ", t= " t ", x= " x " , exact= " (fx t))
    (if (>= i 0)
      (let [h step
            k1 (* h (f x t))
            k2 (* h (f (+ x k1) (+ t h)))
            ex (fx t)]
        (recur (+ x (/ (+ k1 k2) 2)) (+ t h) (dec i)
               (conj result [(double t) (double x) (double ex) (- x ex)])))
      result)))

;;------------------------------------------------------------------------------
;; Runge-Kutta 4th order method
(defn rk4 [f fx x0 t0 step n]
  "Fourth-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0."
  (println "\n\n *** RUNGE-KUTTA 4th ***\n\n")
  (loop [x x0
         t t0
         j n
         result []]
    (println "-- Step= " (- n j) " , t= " t " , x= " x " , exact= " (fx t))
    (if (>= j 0)
      (let [h step
            k1 (* h (f x t))
            k2 (* h (f (+ x (* 0.5 k1)) (+ t (* 0.5 h))))
            k3 (* h (f (+ x (* 0.5 k2)) (+ t (* 0.5 h))))
            k4 (* h (f (+ x k3) (+ t h)))
            ex (fx t)]
        (recur (+ x (/ (+ k1 (* 2 (+ k2 k3)) k4) 6)) (+ t h) (dec j)
               (conj result [(double t) (double x) (double ex) (- x ex)]) ))
      result)))

;;------------------------------------------------------------------------------
;; Runge-Kutta 4th order with error estimate
(defn rk45 [f fx x0 t0 step n]
  "Fourth-order Runge-Kutta method with error estimate.

    USAGE:
        x, err = rk45(f, x0, t)

    INPUT:
        f     - function of x and t equal to dx/dt.  x may be multivalued,
                in which case it should a list or a array.  In this
                case f must return a array with the same dimension
                as x.
        x0    - the initial condition(s).  Specifies the value of x when
                t = t[0].  Can be either a scalar or a list or array
                if a system of equations is being solved.
        t     - list or array of t values to compute solution at.
                t[0] is the the initial condition point, and the difference
                h=t[i+1]-t[i] determines the step size h.

    OUTPUT:
        x     - array containing solution values corresponding to each
                entry in t array.  If a system is being solved, x will be
                an array of arrays.
        err   - array containing estimate of errors at each step.  If
                a system is being solved, err will be an array of arrays.

    NOTES:
        This version is based on the algorithm presented in 'Numerical
        Mathematics and Computing' 6th Edition, by Cheney and Kincaid,
        Brooks-Cole, 2008."
  ;; coeffs used to compute the independent var arg of f
  (def c20 2.500000000000000e-01)
  (def c30 3.750000000000000e-01)
  (def c40 9.230769230769231e-01)
  (def c50 1.000000000000000e+00)
  (def c60 5.000000000000000e-01)
  ;; Coefficients used to compute the dependent variable argument of f
  (def c21 2.500000000000000e-01)
  (def c31 9.375000000000000e-02)
  (def c32 2.812500000000000e-01)
  (def c41 8.793809740555303e-01)
  (def c42 -3.277196176604461e+00)
  (def c43 3.320892125625853e+00)
  (def c51 2.032407407407407e+00)
  (def c52 -8.000000000000000e+00)
  (def c53 7.173489278752436e+00)
  (def c54 -2.058966861598441e-01)
  (def c61 -2.962962962962963e-01)
  (def c62 2.000000000000000e+00)
  (def c63 -1.381676413255361e+00)
  (def c64 4.529727095516569e-01)
  (def c65 -2.750000000000000e-01)
  ;; Coefficients used to compute 4th order RK estimate
  (def a1 1.157407407407407e-01)
  (def a2 0.000000000000000e-00)
  (def a3 5.489278752436647e-01)
  (def a4 5.353313840155945e-01)
  (def a5 -2.000000000000000e-01)
  (def b1 1.185185185185185e-01)
  (def b2 0.000000000000000e-00)
  (def b3 5.189863547758284e-01)
  (def b4 5.061314903420167e-01)
  (def b5 -1.800000000000000e-01)
  (def b6 3.636363636363636e-02)

  (println "\n\n *** RUNGE-KUTTA 45 ***\n\n")
  (loop [x x0
         t t0
         j n
         result []]
    (println "-- Step= " (- n j) " , t= " t " , x= " x " , exact= " (fx t) )
    (if (>= j 0)
      (let [h step
            k1 (* h (f x t))
            x1 (+ x (* c21 k1))
            t1 (+ t (* c20 h))
            k2 (* h (f x1 t1))
            x2 (+ x (+ (* c31 k1) (* c32 k2)))
            t2 (+ t (* c30 h))
            k3 (* h (f x2 t2))
            x3 (+ x (+ (* c41 k1) (* c42 k2) (* c43 k3)))
            t3 (+ t (* c40 h))
            k4 (* h (f x3 t3))
            x4 (+ x (+ (* c51 k1) (* c52 k2) (* c53 k3) (* c54 k4)))
            t4 (+ t h)
            k5 (* h (f x4 t4))
            x5 (+ x (+ (* c61 k1) (* c62 k2) (* c63 k3) (* c64 k4) (* c65 k5)))
            t5 (+ t (* c60 h))
            k6 (* h (f x5 t5))
            x5s (+ x (+ (* b1 k1) (* b3 k3) (* b4 k4) (* b5 k5) (* b6 k6)))
            ex (fx t)
            err (Math/abs (- x5s x))
            ]
        (recur (+ x (+ (* a1 k1) (* a3 k3) (* a4 k4) (* a5 k5)))
               (+ t h)
               (dec j)
               (conj result [t x ex (- x ex)])))
      result)))

;;------------------------------------------------------------------------------
;; Runge-Kutta-Fehlberg method
(defn rkf [f fx a b x0 tol hmax hmin]
  "Runge-Kutta-Fehlberg method to solve x' = f(x,t) with x(t[0]) = x0.

    USAGE:
        t, x = rkf(f, a, b, x0, tol, hmax, hmin)

    INPUT:
        f     - function equal to dx/dt = f(x,t)
        a     - left-hand endpoint of interval (initial condition is here)
        b     - right-hand endpoint of interval
        x0    - initial x value: x0 = x(a)
        tol   - maximum value of local truncation error estimate
        hmax  - maximum step size
        hmin  - minimum step size

    OUTPUT:
        t     - array of independent variable values
        x     - array of corresponding solution function values

    NOTES:
        This function implements 4th-5th order Runge-Kutta-Fehlberg Method
        to solve the initial value problem

           dx
           -- = f(x,t),     x(a) = x0
           dt

        on the interval [a,b].

        Based on pseudocode presented in Numerical Analysis, 6th Edition,
        by Burden and Faires, Brooks-Cole, 1997.
    "
  ;; Coefficients used to compute the independent variable argument of f
  (def a2 2.500000000000000e-01)
  (def a3 3.750000000000000e-01)
  (def a4 9.230769230769231e-01)
  (def a5 1.000000000000000e+00)
  (def a6 5.000000000000000e-01)
  ;; Coefficients used to compute the dependent variable argument of f
  (def b21 2.500000000000000e-01)
  (def b31 9.375000000000000e-02)
  (def b32 2.812500000000000e-01)
  (def b41 8.793809740555303e-01)
  (def b42 -3.277196176604461e+00)
  (def b43 3.320892125625853e+00)
  (def b51 2.032407407407407e+00)
  (def b52 -8.000000000000000e+00)
  (def b53 7.173489278752436e+00)
  (def b54 -2.058966861598441e-01)
  (def b61 -2.962962962962963e-01)
  (def b62 2.000000000000000e+00)
  (def b63 -1.381676413255361e+00)
  (def b64 4.529727095516569e-01)
  (def b65 -2.750000000000000e-01)
  ;; Coefficients used to compute local truncation error estimate.
  (def r1 2.777777777777778e-03 )
  (def r3 -2.994152046783626e-02)
  (def r4 -2.919989367357789e-02)
  (def r5 2.000000000000000e-02)
  (def r6 3.636363636363636e-02)
  ;; Coefficients used to compute 4th order RK estimate
  (def c1 1.157407407407407e-01)
  (def c3 5.489278752436647e-01)
  (def c4 5.353313840155945e-01)
  (def c5 -2.000000000000000e-01)
  ;;
  (println "\n\n**** RUNGE-KUTTA-FEHLBERG ****\n\n")
  (loop [
         x x0
         t a
         j n
         h hmin
         result []]
    (if (< t b)
      (let [
            h (apply min [hmax (- b t)])
            k1 (* h (f x t))
            x1 (+ x (* b21 k1))
            t1 (+ t (* a2 h))
            k2 (* h (f x1 t1))
            x2 (+ x (+ (* b31 k1) (* b32 k2)))
            t2 (+ t (* a3 h))
            k3 (* h (f x2 t2))
            x3 (+ x (+ (* b41 k1) (* b42 k2) (* b43 k3)))
            t3 (+ t (* a4 h))
            k4 (* h (f x3 t3))
            x4 (+ x (+ (* b51 k1) (* b52 k2) (* b53 k3) (* b54 k4)))
            t4 (+ t (* a5 h))
            k5 (* h (f x4 t4))
            x5 (+ x (+ (* b61 k1) (* b62 k2) (* b63 k3) (* b64 k4) (* b65 k5)))
            t5 (+ t (* a6 h))
            k6 (* h (f x5 t5))
            r (/ (Math/abs (+ (* r1 k1) (* r3 k3) (* r4 k4) (* r5 k5) (* r6 k6))) h)
            delta (* 0.84 (Math/pow (/ r tol) 0.25))
            ex (fx t)]
        (println "j= " (- n j) "t= " t ", x= " x  " , ex= " ex ", h= " h
                 "r= " r ", tol= " tol)

        (if (<= r tol)
          (recur (+ x (+ (* c1 k1) (* c3 k3) (* c4 k4) (* c5 k5)))
                 (+ t h)
                 (dec j)
                 (* h (min (max delta 0.1) 4.0))
                 (conj result [t x ex (- x ex)])
                 )
          result)
        )
      result)))

;;------------------------------------------------------------------------------


;; test drive for given f' & f_ex functions
;(println "f= " f " , t0= " t0 " , x0= " x0 " , step= " step " , n= " n)

;; these three lines are for RKF solver
(def tol 1e-6)
(def hmax 0.005)
(def hmin 0.001)



(println "\n--- EULER ---\n")
(euler f exact x0 a b step)

(println "\n--- HEUN METHOD ---\n")
(heun f exact x0 t0 step n)

(println "\n--- ADAMS-BASHFORD ---\n")
(adbash-2 f exact x0 t0 step n)

(println "\n--- ADAMS-BASHFORD-MOULTON PC ---\n")
(abm4-pc f exact x0 t0 step n)

(println "\n--- RK2A ---\n")
(rk2a f exact x0 t0 step n)

(println "\n--- RK2B ---\n")
(rk2b f exact x0 t0 step n)

(println "\n--- RK4 ---\n")
(rk4 f exact x0 t0 step n)

(println "\n--- RK45 ---\n")
(rk45 f exact x0 t0 step n)

(println "\n--- RKF ---\n")
(rkf f exact a b x0 tol hmax hmin)


;(println "\n\n")
;(doseq [q (heun f exact x0 t0 step n)]
;  (println (apply format "t= %.5f f= %.5f exc= %.5f err= %.5f" q)))
;
;(println "\n\n")
;(doseq [q (rk4 f exact x0 t0 step n)]
;  (println (apply format "t= %.5f f= %.5f exc= %.5f err= %.5f" q)))
;
;
;(println "\n\n")
;(doseq [q (rkf f exact a b x0 tol hmax hmin)]
;  (println (apply format "t= %.5f f= %.5f exc= %.5f err= %.5f" q)))

;; exact solution - commented for now
;;(exact b)

(write-to tofile (euler f exact x0 a b step) "\n --- EULER --- \n")
(write-to tofile (heun f exact x0 t0 step n) "\n --- HEUN --- \n")
(write-to tofile (adbash-2 f exact x0 t0 step n) "\n --- ADAMS-BASHFORD --- \n")
(write-to tofile (abm4-pc f exact x0 t0 step n) "\n --- ADAMS-BASHFORD-MOULTON PC --- \n")
(write-to tofile (rk2a f exact x0 t0 step n) "\n --- RK2A --- \n")
(write-to tofile (rk2b f exact x0 t0 step n) "\n --- RK2B --- \n")
(write-to tofile (rk4 f exact x0 t0 step n) "\n --- RK4 --- \n")
(write-to tofile (rk45 f exact x0 t0 step n) "\n --- RK45 --- \n")
(write-to tofile (rkf f exact a b x0 tol hmax hmin) "\n --- RKF --- \n")

