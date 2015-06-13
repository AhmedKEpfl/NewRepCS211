package game;

import java.util.ArrayList;
import java.util.List;

import javax.swing.text.html.HTMLDocument.HTMLReader.BlockAction;

import processing.core.*;

public class Mover extends PApplet {
	
	private PVector location, velocity, acceleration, rotation, friction;
	PApplet app;
	
	private int compteurSmileysHeureux;
	private float score, lastScore;
	private float speed = 1;
	final static float normalForce = 1;
	final static float mu = 0.01f;
	final static float frictionMagnitude = normalForce * mu;
	final static float gravityConstant = 0.1f;
	final static int BLOCK_HEIGHT = 30;
	final static int SCALE_X = 20;
	final static int SCALE_Z = SCALE_X;
	final static float collisionX = BLOCK_HEIGHT * SCALE_X / 2.0f;
	final static float collisionZ = BLOCK_HEIGHT * SCALE_Z / 2.0f;
	final static float BOARD_WIDTH = 2 * collisionX;
	final static float AMORTISSEMENT_COLLISION = 0.5f;
	final static float CYLINDER_BASE_SIZE = 20;
	final static float SHIFT_SCALE_FACTOR = 0.95f, SCALE_FACTOR = 0.6f;
	
	
	List<PVector> listeEnnemis;
	
	PShape arbre;
	
	public void setup(){
		size(TangibleGame.WINDOW_WIDTH, TangibleGame.WINDOW_HEIGHT, P3D);
	}
	
	public Mover(PApplet applet){
		location = new PVector(0 ,0, 0);
		velocity = new PVector(0, 0, 0);
		acceleration = new PVector(0, 0, 0);
		rotation = new PVector(0, 0, 0);
		app = applet;
		
		compteurSmileysHeureux = 0;
		score = 0;
		lastScore = 0;
	    listeEnnemis = new ArrayList<PVector>();
	    listeEnnemis.add(new PVector(200, 0, 200));
	    //barChartColumns = new ArrayList<Integer>();
	    //listeSmileys = new Smiley[numberSmileys];
	    //smileysMechants = new ArrayList<Smiley>();
	    
	    //hs = new HScrollbar(positionScrollX, positionScrollY, 200, 15);
	    
	    //arbre = app.loadShape("C:/Users/Ahmed/Documents/GitHub/assignmentcs211/TGame/src/arbreSimpleFinal.obj");
	    arbre = app.loadShape("arbreSimpleFinal.obj");
	    arbre.rotateZ(PI);
	    
	    
	}
	
	public void update() {
	    friction = velocity.get();
	    friction.mult(-1);
	    friction.normalize();
	    friction.mult(frictionMagnitude);
	    /*if (compteur > COUNTER_BAR_CHART) {
	      barChartColumns.add((int)score);
	      compteur = -1;
	    }
	    compteur++;*/
	    acceleration.x = sin(rotation.z) * gravityConstant;
	    acceleration.z = -sin(rotation.x) * gravityConstant;

	    acceleration.add(friction);

	    velocity.add(acceleration);

	    location.add(velocity);
	    //rectBarChartWidth = RECT_BAR_CHART_WIDTH_INITIAL * hs.getPos();
	  }
	
	public void display() {
	    
	    app.pushMatrix();
	    
	    app.lights();
	    app.stroke(0);
	    app.fill(0);
	    app.background(200, 200, 200);
	    app.text("rotationX: " + degrees(rotation.x) + 
	      "  RotationZ: " + degrees(rotation.z) + 
	      "  Speed: " + speed, 10, 10);
	    app.translate(TangibleGame.WINDOW_WIDTH / 2, TangibleGame.WINDOW_HEIGHT / 2 + BLOCK_HEIGHT);
	    app.scale((TangibleGame.WINDOW_WIDTH / BOARD_WIDTH) * SCALE_FACTOR);

	    //Smileys:
	    /*for(int i = 0; i < numberSmileys; i++){
	      pushMatrix();
	      translate(listeSmileys[i].getX(), listeSmileys[i].getY(), listeSmileys[i].getZ());
	      scale(listeSmileys[i].getScale());
	      listeSmileys[i].lookAt(location.x, location.y, location.z);
	      rotateY(listeSmileys[i].getRotationY());
	      rotateX(listeSmileys[i].getRotationX());
	      if(listeSmileys[i].getHeureux()){
	        shape(smileyHeureux);
	      } else {
	        shape(smileyTriste);
	      }
	      popMatrix();
	    }*/
	    app.rotateX(rotation.x - PI / 5);
	    app.rotateY(rotation.y);
	    app.rotateZ(rotation.z);
	    app.pushMatrix();
	    if(degrees(rotation.x) > 36){
	    	app.pushMatrix();
	    	app.translate(location.x, -BLOCK_HEIGHT, location.z);
	    	app.stroke(0);
	    	app.fill(0);
	    	app.sphere(BLOCK_HEIGHT / 2.0f);
	    	
	    	app.popMatrix();
	    	for (int i = 0; i < listeEnnemis.size (); i++) {
		  	      displayEnnemi(listeEnnemis.get(i));
		  	    }
	    }
	    app.scale(SCALE_X, 1, SCALE_Z);
	    app.stroke(50, 50, 50, 50);
	    app.fill(50, 50, 50, 50);
	    app.box(BLOCK_HEIGHT);
	    app.popMatrix();
	    for (int i = 0; i < listeEnnemis.size (); i++) {
	      displayEnnemi(listeEnnemis.get(i));
	    }

	    app.translate(location.x, -BLOCK_HEIGHT, location.z);
	    app.stroke(100, 100, 100);
	    app.fill(200, 200, 200);
	    app.sphere(BLOCK_HEIGHT / 2.0f);
	    app.popMatrix();
	    /*app.drawDataVisualization();
	    image(dataVisualization, 0, height - dataVisualization.height);
	    hs.update();
	    hs.display();*/
	    
	}
	
	public void displayShift() {
	    app.lights();
	    
	    app.translate(TangibleGame.WINDOW_WIDTH / 2f, TangibleGame.WINDOW_HEIGHT / 2f, 0);
	    app.rotateX(-PI / 2);
	    app.scale(SHIFT_SCALE_FACTOR * TangibleGame.WINDOW_WIDTH / BOARD_WIDTH);
	    app.pushMatrix();
	    app.scale(SCALE_X, 1, SCALE_Z);

	    app.stroke(50, 50, 50);
	    app.fill(50, 50, 50);
	    app.box(BLOCK_HEIGHT);
	    app.popMatrix();
	    for (int i = 0; i < listeEnnemis.size(); i++) {
	      displayEnnemi(listeEnnemis.get(i));
	    }
	    app.translate(location.x, -BLOCK_HEIGHT, location.z);
	    app.stroke(100, 100, 100);
	    app.fill(200, 200, 200);
	    app.sphere(BLOCK_HEIGHT / 2.0f);
	  }
	
	public void checkEdges(){
		//Check x-colision:
		if (location.x >= collisionX){
			if(velocity.x > 0.1){
				changeScoreBorder();
			}
			location.x = collisionX;
			velocity.x = -velocity.x * AMORTISSEMENT_COLLISION;
		} else if(location.x <= -collisionX){
			if(velocity.x < -0.1){
				changeScoreBorder();
			}
			location.x = -collisionX;
			velocity.x = -velocity.x * AMORTISSEMENT_COLLISION;
		}
		
		//Check z-collision:
		if (location.z >= collisionZ){
			if(velocity.z > 0.1){
				changeScoreBorder();
			}
			location.z = collisionZ;
			velocity.z = -velocity.z * AMORTISSEMENT_COLLISION;
		} else if(location.z <= -collisionZ){
			if(velocity.z < -0.1){
				changeScoreBorder();
			}
			location.z = -collisionZ;
			velocity.z = -velocity.z * AMORTISSEMENT_COLLISION;
		}
	}
	
	public void changeScoreBorder(){
		score -= velocity.mag();
		lastScore = velocity.mag();
		if(score < 0){
			score = 0;
		}
	}
	
	public void checkCylinderCollision() {
	    int cylinderCollision = -1;
	    PVector n = new PVector(0, 0, 0);
	    PVector nCopy = new PVector(0, 0, 0);

	    for (int i = 0; i < listeEnnemis.size (); i++) {
	      if (CYLINDER_BASE_SIZE + BLOCK_HEIGHT / 2.0 >= listeEnnemis.get(i).dist(location)) {
	        /*if(compteurSmileysHeureux < NUMBER_SMILEYS){
	          listeSmileys[compteurSmileysHeureux].setHeureux();
	          compteurSmileysHeureux++;
	        }*/
	        
	        score += velocity.mag();
	        lastScore = velocity.mag();
	        cylinderCollision = i;
	        
	        n.x = (location.x - (listeEnnemis.get(i)).x);
	        n.y = (location.y - (listeEnnemis.get(i)).y);
	        n.z = (location.z - (listeEnnemis.get(i)).z);
	        n.normalize();
	        nCopy = n.get();

	        n.mult(2* (velocity.dot(n)));
	        velocity.sub(n);
	        nCopy.mult(CYLINDER_BASE_SIZE + BLOCK_HEIGHT / 2.0f);
	        location.set(nCopy);
	        location.add(listeEnnemis.get(i));
	        listeEnnemis.remove(i);
	        //n * r +  + listeCylindre.get(i)
	      }
	    }
	  }
	
	public void setRotationX(float newValue){
		if(newValue > PI / 3){
			rotation.x = PI / 3;
		} else if(newValue < -PI / 3){
			rotation.x = -PI / 3;
		} else{
			rotation.x = newValue;
		}
	}
	
	public void setRotationZ(float newValue){
		if(newValue > PI / 3){
			rotation.z = PI / 3;
		} else if(newValue < -PI / 3){
			rotation.z = -PI / 3;
		} else{
			rotation.z = newValue;
		}
	}
	
	public void addRotationX(float addValue){
		rotation.x += addValue;
		if(rotation.x > PI / 3){
			rotation.x = PI / 3;
		} else if(rotation.x < -PI / 3){
			rotation.x = -PI / 3;
		}
	}
	
	public void addRotationZ(float addValue){
		rotation.z += addValue;
		if(rotation.z > PI / 3){
			rotation.z = PI / 3;
		} else if(rotation.z < -PI / 3){
			rotation.z = -PI / 3;
		}
	}
	
	public float getRotationX(){
		return rotation.x;
	}
	
	public float getRotationZ(){
		return rotation.z;
	}
	
	public float getSpeed(){
		return speed;
	}
	
	public void displayEnnemi(PVector ennemi){
		app.pushMatrix();
		app.translate(ennemi.x,  - 2 * BLOCK_HEIGHT, ennemi.z);
		app.shape(arbre);
		
		app.popMatrix();
	}
	
	public void addEnnemi(PVector ennemi){
		listeEnnemis.add(new PVector(ennemi.x, ennemi.y, ennemi.z));
	}
}
