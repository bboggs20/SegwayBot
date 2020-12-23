#include <EncoderStepCounter.h>

#include <SparkFun_TB6612.h>


#define AIN1 5
#define AIN2 4
#define PWMA 6
#define STBY 7
#define ENCA 8
#define ENCB 9

#define EINTA digitalPinToInterrupt(ENCA)
#define EINTB digitalPinToInterrupt(ENCB)


#define dt 20 // velocity measurement dt in ms

EncoderStepCounter encoder(ENCA, ENCB);
Motor q1 = Motor(AIN1, AIN2, PWMA, 1, STBY);

long position, pos0;
double velocity;

unsigned long t0, t;

void setup() 
{

    position = 0; pos0 = 0;
    velocity = 0;
    t0 = millis(); 
    t = millis();
    
    Serial.begin(9600);
    encoder.begin();
    
    attachInterrupt(EINTA, interrupt, CHANGE);
    attachInterrupt(EINTB, interrupt, CHANGE);
    
    q1.drive(255);
    //move_to(q1,encoder,&position,1000);
    
}

void interrupt()
{
  encoder.tick();
}

void loop() 
{
  updateEncState(encoder, &position);
  updateVelocity();
  
  Serial.print("Pos: ");
  Serial.println(position);
  Serial.print("Vel: ");
  Serial.println(velocity);
}


void move_to(Motor m, EncoderStepCounter enc, long* p, long pos_d)
{
  int kp = 100, kd = 10;
  int pwmh[2] = {0,0};
  int v = 0, pwm = 0;
  unsigned long t0 = millis(), t = 0;
  while (*p != pos_d)
  {
    t = millis();
    v = (pwmh[1]-pwmh[0])/(t-t0);
    t0 = t;
    pwm = constrain(-kp*(*p-pos_d)-kd*v, -255, 255);
    m.drive(pwm);
    pwmh[1] = pwmh[0];
    pwmh[0] = pwm;
    updateEncState(enc, p);
  }
}

void updateVelocity()
{
  t = millis();
  if (t-t0 >= dt)
  {
    velocity = (1.0*(position-pos0))/(t-t0);
    t0 = t;
    pos0 = position;
  }
}
void updateEncState(EncoderStepCounter enc, long* p)
{
  char pos = enc.getPosition();
  if (pos != 0) {
    *p += pos;
    encoder.reset();
    //Serial.println(position);
  }
}
