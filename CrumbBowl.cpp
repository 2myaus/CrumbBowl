#include <mutex>
#include <vector>
#include <array>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <execution>
#include <chrono>
#include <cstdlib>

class Vector2{
  public:
  double x;
  double y;

  Vector2(double argX, double argY): x(argX), y(argY) {}
  Vector2(): x(0), y(0) {}

  Vector2 operator+(const Vector2 &o){
    return Vector2(x+o.x, y+o.y);
  }
  Vector2 operator-(const Vector2 &o){
    return Vector2(x-o.x, y-o.y);
  }

  Vector2 operator*(const double &o){
    return Vector2(x*o, y*o);
  }
  Vector2 operator/(const double &o){
    return Vector2(x/o, y/o);
  }

  void operator+=(const Vector2 &o){
    this->x += o.x;
    this->y += o.y;
  }

  void operator-=(const Vector2 &o){
    this->x -= o.x;
    this->y -= o.y;
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
  double getMinPos() { return minPos; }
  double getMaxPos() { return maxPos; }

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
      creating.grid[gridY][gridX].totalMass += particle.mass;
    });

    return creating;
  }

  static ParticleGroupGrid fromHigherDetailGrid(ParticleGroupGrid &argOtherGrid){
    assert(argOtherGrid.detail > 0);

    ParticleGroupGrid creating;
    creating.detail = argOtherGrid.getDetail()-1;

    const unsigned int gridWidth = std::pow(2, creating.detail);

    creating.grid.resize(gridWidth);
    std::for_each(std::execution::seq, creating.grid.begin(), creating.grid.end(), [&](std::vector<ParticleGroup> &gridCol){
      gridCol.resize(gridWidth);
    });

    for(unsigned int row = 0; row < gridWidth; row++){
      for(unsigned int col = 0; col < gridWidth; col++){
        creating.grid.at(row).at(col).subGroups[0] = &argOtherGrid.grid[row*2][col*2];
        creating.grid[row][col].subGroups[1] = &argOtherGrid.grid[row*2][col*2+1];
        creating.grid[row][col].subGroups[2] = &argOtherGrid.grid[row*2+1][col*2];
        creating.grid[row][col].subGroups[3] = &argOtherGrid.grid[row*2+1][col*2+1];

        creating.grid[row][col].totalMass =
          creating.grid[row][col].subGroups[0]->totalMass +
          creating.grid[row][col].subGroups[1]->totalMass +
          creating.grid[row][col].subGroups[2]->totalMass +
          creating.grid[row][col].subGroups[3]->totalMass;
      }
    }

    creating.minPos = argOtherGrid.minPos;
    creating.maxPos = argOtherGrid.maxPos;
    
    return creating;
  }
};

class ParticleGroupGridCollection{
  std::vector<ParticleGroupGrid> detailLevels;
  int maxDetailLevel;

  public:
  ParticleGroupGrid *getDetailLevel(int detailLevel){
    return &(detailLevels.at(detailLevel));
  }

  const int getMaxDetail(){ return maxDetailLevel; }

  static ParticleGroupGridCollection fromParticles(std::vector<Particle> &argParticles, int argMaxDetail){
    ParticleGroupGridCollection creating;
    creating.maxDetailLevel = argMaxDetail;

    creating.detailLevels.resize(argMaxDetail+1);
    creating.detailLevels[argMaxDetail] = ParticleGroupGrid::fromParticles(argParticles, argMaxDetail);
    for(int detail = argMaxDetail - 1; detail >= 0; detail--){
      creating.detailLevels[detail] = ParticleGroupGrid::fromHigherDetailGrid(creating.detailLevels[detail+1]);
    }

    return creating;
  }
};

Vector2 getAttraction(Vector2 dif, double m1, double m2){ // dv/dt
  return Vector2(); //TODO
}

void simulateTheseParticles(std::vector<Particle*> &theseParticles, std::vector<Particle*> &consideringThese, Vector2 accelOnAll, double deltaTime){
  std::for_each(std::execution::par, theseParticles.begin(), theseParticles.end(), [&](Particle *simulatedParticle){
    std::for_each(std::execution::par, theseParticles.begin(), theseParticles.end(), [&](Particle *consideredParticle){
      if(simulatedParticle == consideredParticle){ return; }
      Vector2 accel = getAttraction(consideredParticle->position-simulatedParticle->position, simulatedParticle->mass, consideredParticle->mass) + accelOnAll;
      simulatedParticle->position += (simulatedParticle->velocity * deltaTime) + (accel*0.5*deltaTime*deltaTime);
      simulatedParticle->velocity += accel*deltaTime;
    });
    std::for_each(std::execution::par, consideringThese.begin(), consideringThese.end(), [&](Particle *consideredParticle){
      Vector2 accel = getAttraction(consideredParticle->position-simulatedParticle->position, simulatedParticle->mass, consideredParticle->mass) + accelOnAll;
      simulatedParticle->position += (simulatedParticle->velocity * deltaTime) + (accel*0.5*deltaTime*deltaTime);
      simulatedParticle->velocity += accel*deltaTime;
    });
  });
}

void simulateParticles(ParticleGroupGridCollection &particleCollection, double deltaTime){
  std::vector<std::vector<Vector2>> velocityPrecalcs;

  {
    unsigned int highestSquare = std::pow(2, particleCollection.getMaxDetail());
    highestSquare *= highestSquare;
    velocityPrecalcs.resize(particleCollection.getMaxDetail()+1);
    std::for_each(velocityPrecalcs.begin(), velocityPrecalcs.end(), [&](auto &velocityVec){
      velocityVec.resize(highestSquare);
    });
  }

  for(int detailLevel = 0; detailLevel <= particleCollection.getMaxDetail(); detailLevel++){
    unsigned int gridWidth = std::pow(2, detailLevel);
    ParticleGroupGrid *currentGrid = particleCollection.getDetailLevel(detailLevel);

    double currentBoundsSize = currentGrid->getMaxPos() - currentGrid->getMinPos();

    for(int row = 0; row < gridWidth; row++){
      for(int col = 0; col < gridWidth; col++){
        int localRowMid = std::floor(row / 2) * 2;
        int localColMid = std::floor(col / 2) * 2;
        //Apply forces for squares only within [localRowMid-2, localRowMid + 4) (and same for col)
        
        for(int row2 = localRowMid-2; row2 < localRowMid+4; row2++){
          for(int col2 = localColMid-2; col2 < localColMid+4; col2++){
            if(row2 < 0 || col2 < 0 || row2 >= gridWidth || col2 >= gridWidth){ continue; }

            if((
              row2+1 == row ||
              row2 == row ||
              row == row+1
              ) && (
              col2+1 == col ||
              col2 == col ||
              col2 == col+1
              ))
            {
              //The two squares are in the 3x3 area around each other
              continue;
            }

            Vector2 squaresDist = (Vector2(col2, row2) - Vector2(col, row)) * currentBoundsSize;

            velocityPrecalcs[detailLevel][row*gridWidth+col] += getAttraction(
                squaresDist,
                currentGrid->getGroupAt(row, col)->totalMass,
                currentGrid->getGroupAt(row2, col2)->totalMass
            );

            //TODO: add the velocityPrecalcs to the higher-detail groups below, unless on highest-detail grid
            if(detailLevel == particleCollection.getMaxDetail()){ continue; }
            velocityPrecalcs.at(detailLevel+1).at(row*gridWidth*4+col*2) += velocityPrecalcs.at(detailLevel).at(row*gridWidth+col);
            velocityPrecalcs.at(detailLevel+1).at((row*2+1)*gridWidth*2+col*2) += velocityPrecalcs[detailLevel][row*gridWidth+col];
            velocityPrecalcs.at(detailLevel+1).at((row*2+1)*gridWidth*2+(col*2)+1) += velocityPrecalcs[detailLevel][row*gridWidth+col];
            velocityPrecalcs.at(detailLevel+1).at(row*gridWidth*4+(col*2)+1) += velocityPrecalcs[detailLevel][row*gridWidth+col];
          }
        }

        //Apply particle-to-particle forces for the 3x3 grid around this point, if we are at the highest detail level
        if(detailLevel != particleCollection.getMaxDetail()){ continue; }
        
        std::vector<Particle*> localParticles;
        localParticles.reserve(10000); //TODO: Calculate/reserve this dynamically based on stats

        for(int row2 = row-1; row2 <= row+1; row2++){
          for(int col2 = col-1; col2 <= col+1; col2++){
            if(row2 == row && col2 == col){ continue; }
            if(row2 < 0 || col2 < 0 || row2 >= gridWidth || col2 >= gridWidth){ continue; }

            auto currentGridGroup = currentGrid->getGroupAt(row2, col2);
            localParticles.insert(localParticles.end(), currentGridGroup->particles.begin(), currentGridGroup->particles.end());
          }
        }

        simulateTheseParticles(currentGrid->getGroupAt(row, col)->particles, localParticles, velocityPrecalcs[detailLevel][row*gridWidth+col], deltaTime);
      }
    }
  }
}

int main(){
  unsigned int num = 200000;
  std::vector<Particle> particles;
  particles.reserve(num);

  for(unsigned int i = 0; i < num; i++){
    particles.push_back(Particle::random());
  }


  auto begin = std::chrono::steady_clock::now();
  ParticleGroupGridCollection collection = ParticleGroupGridCollection::fromParticles(particles, 7);
  auto end1 = std::chrono::steady_clock::now();
  simulateParticles(collection, 0.1);
  auto end2 = std::chrono::steady_clock::now();


  printf("Placement on grid took %ld ns!", std::chrono::duration_cast<std::chrono::nanoseconds> (end1 - begin).count());
  printf("Simulation took %ld ns!", std::chrono::duration_cast<std::chrono::nanoseconds> (end2 - end1).count());
}
