# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

---

## Notes

### Radar not tracking

This seems like a bug... without lidar, it is pretty hopeless on dataset 2.
It does not too badly on dataset 1...
And I guess on 2 it may be more to do with numerical failure.

### Ambiguity in v / phi

The speed can go negative, which represents the same state as `-v` / `phi + pi`.

### Numerical issues due to large yaw rates

If the yaw rate gets large, the sigma points generate essentially random yaws --- basically spinning a wheel and seeing where it stops. This causes fluctuations in yaw, which means the error signal on the yaw rate contains no information, and it is underdamped.

I think the Gaussian assumption is only valid when both the yaw rate and the yaw rate variance are small --- if the target is spinning rapidly, we can't tell whether it's rotating at xHz or 2xHz. So... we could normalise the yaw rate as well as the yaw. That would prevent it from blowing up.

Do we also need to do something about its variance, to prevent the variance from blowing up? In reality, if we have variance of ~2pi, we have no idea what the yaw rate is. We could just clamp it, but...

Actually, I think the problem is that we get a negative variance term in P, because lambda is too negative.

### Output does not tell you whether it's laser or radar NIS

So you can't really tell whether you have 2 or 3 degrees of freedom. Maybe need a script to augment, or a flag to change the output format to include both NIS values.

### Need to tune parameters

Just 3 main ones:

- process noise in dv
- process noise in dyaw
- lambda

A few hidden parameters:
- initial noise values
- div by zero tolerance
- if we do more aggressive normalisation, there may be more

Setting a low value acceleration noise value seems to be required for dataset 2. A lower yaw noise also helps on dataset 2, but setting the yaw noise too low increases RMSE on dataset 1.

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) as much as possible.

## Generating Additional Data

This is optional!

If you'd like to generate your own radar and lidar data, see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.

## Project Instructions and Rubric

This information is only accessible by people who are already enrolled in Term 2
of CarND. If you are enrolled, see [the project page](https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/c3eb3583-17b2-4d83-abf7-d852ae1b9fff/concepts/4d0420af-0527-4c9f-a5cd-56ee0fe4f09e)
for instructions and the project rubric.
