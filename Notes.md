# Overall comments
The project compiles and runs without issues for the graphic-less path. Had to install many packages for OpenGL (not specified by the readme) but we succeeded in running it. The visual results are really nice. Moreover there is a good usage of git submodules, often there are useful comments and the CMake provides a release and debug profile with appropriate flags (a pity that it is not mentioned in the readme).

# Technical comments
* In modern C++ instead of
```cpp
Vector<DIM> position, speed, acceleration;
double mass(static_cast<double>((std::rand() % 10000) * 1.0e7 + 1.0e7));
position = Vector<DIM>({static_cast<double>(std::rand() % 4000) - 2000.0, static_cast<double>(std::rand() % 4000) - 2000.0, static_cast<double>(std::rand() % 4000) - 2000.0});
speed = Vector<DIM>({0.0, 0.0, 0.0});
acceleration = Vector<DIM>({0.0, 0.0, 0.0});
```
it would be more idiomatic something like
```cpp
Vector<DIM> position, speed, acceleration;
const auto mass = static_cast<double>((std::rand() % 10000) * 1.0e7 + 1.0e7);
const Vector<DIM> position = {{static_cast<double>(std::rand() % 4000) - 2000.0, static_cast<double>(std::rand() % 4000) - 2000.0, static_cast<double>(std::rand() % 4000) - 2000.0}};
const Vector<DIM> speed = {{0.0, 0.0, 0.0}};
const Vector<DIM> acceleration = {{0.0, 0.0, 0.0}};
```
moreover, you should rely on `<random>` instead of `std::rand`.
* Be careful in using `using namespace std;` It is relatively safe to do it only inside the body of the `main()`, otherwise it is dangerous since you are bringing to the local
scope all names of the included headers of the standard library! If you wish to just bring a 
name into the current scope, for the sake of simplicity, you can do, for nstace `using std::vector`.
* The parameters defined in the global scope in `main.cpp` should be put in a separate configuration file.
* In `json_parse.cpp` sometimes you pass strings by copy instead of const reference. Use `const &` for potentially "large" objects.
* I do not see why `global_particles` and `global_relocated_particles` are global variables. If there is not a very specific reason for them to be global this is bad practice. Avoid global variables as far as possible!
* `sumUpTo` could be replaced using the explicit formula for the sum for the first n integers.
* In `constexpr friend Vector<dim> operator+(const Vector<dim> lhs, const Vector<dim>& rhs)` you pass lhs by copy, but you do not exploit the fact that it is a copy and you copy it again inside the body of the function. By the way, why `constexpr`? It is true that nowadays constexpr is popping up everywhere, but here I do not see how the compiler can take advantage of it since you cannot have a constexpr `std::vector` (at least it does not make sense to me). Anyway, it is not an error, if the function cannot be executed at compile time `constexpr` is equivalent to `inline`. Similar comments apply to `__create_subcluster_at`, `_get_cluster_subscript`
* Don't take the habit of starting a name with one or two underscores. Normally starting with underscore is used by system constants, and you do not want the risk of overriding a system constant by mistake. Underscores are fine in the middle or at the end.
* In `pow(distance, 3)` you are not using the `std` version of the function and may be inefficient when the exponent is a small integer, better `distance*distance*distance` (may be untrue in latest version of the compiler, but only if you use the `std::` version).
* `Particle` constructors takes copies instead of const references, moreover `bodies` is unused. Take the habit of using  `const &` instead of passing by value, particularly for object that may be big.
  With the possible exception of int and double.
* You compute `std::array<unsigned int, dim> dim_mods` at each call of `ParticleCluster<dim>::_create_subclusters_one_level`, this should have been a `static` variable, or even a `constexpr static`
* Why do you use `std::map<unsigned int, ParticleCluster<dim>> children_clusters;` instead of `std::vector`? It seems to me that the key of the map is a range (it is not sparse) so a vector would be more efficient


# Minor issues
* In `json_parse.cpp` an `unordered_map` should perform better than a `map`
* Avoid folder (and file) names with spaces!
* Should remove the unused file `main_old.cpp`, moreover never put binaries in repos (see `main_old`)
* Code that is not used should be removed from the git repo (don't worry, git has good memory, you can still recover them if needed).
* There are some C-style cast for `double`s. It is not an error, but in C++ it is better to use the C++-style casts because they are safer.
* I would argue that public inheritance is not the best design pattern for `ParticleCluster`, probably composition makes more sense.
* `local_max_boundary = Particle<dim>::max_boundary + Particle<dim>::max_boundary * .01;` is such a rough approximation needed? Is it not enough to use a small approximation like 1e-8 instead of 0.01?


