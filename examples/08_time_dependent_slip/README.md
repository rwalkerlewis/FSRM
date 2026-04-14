# Example 08: Time-dependent Prescribed Slip

## Physics

Applies a prescribed strike-slip displacement on a vertical fault plane with a
linear time ramp. The slip ramps from zero to 1 m over the interval [0.1, 0.6] seconds,
then remains constant at 1 m for the rest of the simulation.

This exercises the time-dependent prescribed slip feature:
- slip(t) = 0 for t < onset_time
- slip(t) = slip_max * (t - onset_time) / rise_time for onset_time <= t < onset_time + rise_time
- slip(t) = slip_max for t >= onset_time + rise_time

## Config

`config/examples/time_dependent_slip.config`

## Fault Parameters

| Parameter          | Value   | Description                     |
|--------------------|---------|---------------------------------|
| mode               | prescribed_slip | Kinematic fault |
| strike             | 0 deg   | N-S fault plane                 |
| dip                | 90 deg  | Vertical fault                  |
| slip_strike        | 1.0 m   | Maximum right-lateral slip      |
| slip_onset_time    | 0.1 s   | Time when slip begins ramping   |
| slip_rise_time     | 0.5 s   | Duration of linear ramp         |

## Running

```bash
./run.sh
```

Or manually:

```bash
cd build
./fsrm -c ../config/examples/time_dependent_slip.config
```

## Expected Output

- `output/solution.h5`: Time series of displacement field
- At t = 0.1: slip = 0
- At t = 0.35: slip = 0.5 m (50% of max)
- At t = 0.6: slip = 1.0 m (full slip)
- At t = 1.0: slip = 1.0 m (unchanged)

## Verified By

- `Integration.TimeDependentSlip` test
