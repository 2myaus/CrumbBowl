#include <future>
#include <mutex>
#include <vector>
#include <array>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <memory>
#include <execution>
#include <chrono>
#include <cstdlib>

class Vector2{
  public:
  double x;
  double y;

  Vector2(double argX, double argY): x(argX), y(argY) {}
  Vector2(): x(0), y(0) {}

  Vector2 operator+(Vector2 &o){
    return Vector2(x+o.x, y+o.y);
  }
  Vector2 operator-(Vector2 &o){
    return Vector2(x-o.x, y-o.y);
  }

  Vector2 operator*(double &o){
    return Vector2(x*o, y*o);
  }
  Vector2 operator/(double &o){
    return Vector2(x/o, y/o);
  }

};

class Particle{
  public:
  double mass;
  Vector2 position;
  Vector2 velocity;

  Particle(Vector2 argPosition, Vector2 argVelocity): position(argPosition), velocity(argVelocity){}

  static Particle random(){
    return Particle(Vector2((double)rand() / RAND_MAX, (double)rand() / RAND_MAX), Vector2((double)rand() / RAND_MAX, (double)rand() / RAND_MAX));
  }
};

class ParticleGroup{
  public:
  std::vector<Particle*> particles; // When used in the highest-detail grid, lists particles that fall within this group
  std::array<ParticleGroup*, 4> subGroups; // When used in a lower-detail grid, refers to corresponding ParticleGroups on the higher detail grid

  double totalMass;
  Vector2 averagePosition;
  Vector2 averageVelocity;

  ParticleGroup(): totalMass(0), averagePosition(), averageVelocity(){}
};

class ParticleGroupGrid{
  double maxPos;
  double minPos;
  unsigned int detail;
  std::vector<std::vector<ParticleGroup>> grid;

  public:
  unsigned int getDetail() { return detail; }
  ParticleGroup *getGroupAt(unsigned int row, unsigned int col){
    assert(row < grid.size() && col < grid.at(row).size());

    return &(grid.at(row).at(col));
  }

  ParticleGroupGrid(): maxPos(0), minPos(0){};

  static ParticleGroupGrid fromParticles(std::vector<Particle> &argParticles, unsigned int argDetail){
    ParticleGroupGrid creating;
    creating.detail = argDetail;

    unsigned int gridWidth = std::pow(2, argDetail);

    std::for_each(argParticles.begin(), argParticles.end(), [&](Particle &particle){
      if(particle.position.x < creating.minPos){ creating.minPos = particle.position.x; }
      else if(particle.position.x > creating.maxPos) { creating.maxPos = particle.position.x; }

      if(particle.position.y < creating.minPos){ creating.minPos = particle.position.y; }
      else if(particle.position.y > creating.maxPos) { creating.maxPos = particle.position.y; }
    });

    double size = (creating.maxPos-creating.minPos);

    auto mutexGrid = std::vector<std::mutex>(gridWidth*gridWidth);

    creating.grid.resize(gridWidth);
    std::for_each(std::execution::seq, creating.grid.begin(), creating.grid.end(), [&](std::vector<ParticleGroup> &gridCol){
      gridCol.resize(gridWidth);
    });

    std::for_each(std::execution::par, argParticles.begin(), argParticles.end(), [&](Particle &particle){
      double xFromMin = particle.position.x - creating.minPos;
      double yFromMin = particle.position.y - creating.minPos;

      if(xFromMin < 0) xFromMin = 0;
      if(yFromMin < 0) yFromMin = 0;

      unsigned int gridX = std::floor(xFromMin / size);
      unsigned int gridY = std::floor(yFromMin / size);

      if(gridX >= gridWidth) { gridX = gridWidth - 1; }
      if(gridY >= gridWidth) { gridY = gridWidth - 1; }

      std::lock_guard<std::mutex> lock(mutexGrid.at(gridY*gridWidth+gridX));
      creating.grid.at(gridY).at(gridX).particles.push_back(&particle);
    });

    return creating;
  }
};

int main(){
  unsigned int num = 20000000;
  std::vector<Particle> particles;
  particles.reserve(num);

  for(unsigned int i = 0; i < num; i++){
    particles.push_back(Particle::random());
  }


  auto begin = std::chrono::steady_clock::now();
  ParticleGroupGrid::fromParticles(particles, 5);
  auto end = std::chrono::steady_clock::now();

  printf("Took %ld ns!", std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count());
}
