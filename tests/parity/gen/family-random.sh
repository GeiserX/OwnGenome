# shellcheck shell=bash
# The fixed-seed random stream used when these fixtures were frozen.
# All intermediate integers fit exactly in an awk double.
# shellcheck disable=SC2034 # The family generators pass this program to awk.
family_random=$(cat <<'AWK'
function family_seed(seed,     i) {
  family_rng[0] = seed
  for (i = 1; i < 31; i++) family_rng[i] = (16807 * family_rng[i - 1]) % 2147483647
  for (i = 31; i < 34; i++) family_rng[i % 31] = family_rng[(i - 31) % 31]
  family_index = 34
  while (family_index < 344) family_rand()
}
function family_rand(     i, value) {
  i = family_index++
  value = (family_rng[(i - 31) % 31] + family_rng[(i - 3) % 31]) % 4294967296
  family_rng[i % 31] = value
  return int(value / 2) / 2147483648
}
AWK
)
