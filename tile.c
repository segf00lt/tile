#include "basic.h"
#include <raylib.h>
#include <raymath.h>
#include <math.h>
#define STB_SPRINTF_IMPLEMENTATION
#include "stb_sprintf.h"

const float ENTITY_MOVE_FORCE = 7e5;
const float FRICTION = 18.20;

const float TILE_SIZE = 32.0;
const float HALF_TILE_SIZE = TILE_SIZE * 0.5;
#define NEGATIVE_ZERO (u32)(0x80000000)
#define CHUNK_WIDTH 32
#define CHUNK_HEIGHT 32
u8 CHUNK[CHUNK_HEIGHT * CHUNK_WIDTH] = {
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

typedef struct Entity Entity;

INLINE bool checkCollisionPointRecNoBorder(Vector2 point, Rectangle rec);
INLINE bool lineSegIntersect(Vector2 p1, Vector2 p2, Vector2 p3, Vector2 p4, Vector2 *p);
INLINE void canonicalize(Vector2 p, u64 *x, u64 *y);
INLINE void entityGetInput(Entity *ent);
INLINE void entitySimulate(Entity *ent, float timestep);
INLINE void entityDraw(Entity *ent);

struct Entity {
	Vector2 force;
	Vector2 accel;
	Vector2 vel;
	Vector2 pos;
	float mass;
	float invmass; // 1/mass
	float widthH; // half width
	float heightH; // half height
};

INLINE bool lineSegIntersect(Vector2 p1, Vector2 p2, Vector2 p3, Vector2 p4, Vector2 *p) {
	float t, u;
	float numerator, denominator;

	float p12x = p1.x - p2.x;
	float p12y = p1.y - p2.y;
	float p34x = p3.x - p4.x;
	float p34y = p3.y - p4.y;
	float p13x = p1.x - p3.x;
	float p13y = p1.y - p3.y;

	denominator = p12x*p34y - p12y*p34x;
	denominator = 1/denominator;

	numerator = p13x * p34y - p13y * p34x;

	t = numerator * denominator;

	numerator = p13x * p12y - p13y * p12x;

	u = numerator * denominator;

	*p = (Vector2){ p1.x + t*(p2.x - p1.x), p1.y + t*(p2.y - p1.y) };

	return (0.0 <= t && t <= 1.0) && (0.0 <= u && u <= 1.0);
}

INLINE bool checkCollisionPointRecNoBorder(Vector2 point, Rectangle rec) {
	bool collision = false;

	if((point.x>rec.x)&&(point.x<(rec.x+rec.width))&&(point.y>rec.y)&&(point.y<(rec.y+rec.height)))
		collision = true;

	return collision;
}

INLINE void canonicalize(Vector2 p, u64 *x, u64 *y) {
	*x = Clamp((float)(u64)(p.x) / TILE_SIZE, 0, CHUNK_WIDTH);
	*y = Clamp((float)(u64)(p.y) / TILE_SIZE, 0, CHUNK_HEIGHT);
}

INLINE void entityGetInput(Entity *ent) {
	ent->force = (Vector2){ 0, 0 };

	if(IsKeyDown(KEY_UP)) {
		ent->force.y = -ENTITY_MOVE_FORCE;
	}

	if(IsKeyDown(KEY_LEFT)) {
		ent->force.x = -ENTITY_MOVE_FORCE;
	}

	if(IsKeyDown(KEY_RIGHT)) {
		ent->force.x = ENTITY_MOVE_FORCE;
	}

	if(IsKeyDown(KEY_DOWN)) {
		ent->force.y = ENTITY_MOVE_FORCE;
	}
}

INLINE bool rayAABBIntersect(Vector2 a, Vector2 invdir, Vector2 min, Vector2 max, float *toi) {
	float t1, t2, tmin, tmax;
    t1 = (min.x - a.x) * invdir.x;
    t2 = (max.x - a.x) * invdir.x;
    tmin = fminf(t1, t2);
    tmax = fmaxf(t1, t2);
	t1 = (min.y - a.y) * invdir.y;
	t2 = (max.y - a.y) * invdir.y;
	tmin = fmaxf(tmin, fminf(fminf(t1, t2), tmax));
	tmax = fminf(tmax, fmaxf(fmaxf(t1, t2), tmin));
	*toi = (tmin < 0) ? tmax : tmin;
	return tmax > fmaxf(tmin, 0.0);
}

INLINE void entitySimulate(Entity *ent, float timestep) {
	/*
	 * a = f/m
	 * v = a*t + v0
	 * s = 0.5*a*t^2 + v*t + s0
	 */

	ent->accel = Vector2Scale(ent->force, ent->invmass * timestep);
	ent->vel = Vector2Subtract(Vector2Add(ent->vel, ent->accel), Vector2Scale(ent->vel, timestep*FRICTION));
	Vector2 offset = Vector2Add(Vector2Scale(ent->vel,timestep),Vector2Scale(ent->accel,HALF(timestep)));
	Vector2 newpos = Vector2Add(ent->pos, offset);

	u64 canonx1, canony1;
	u64 canonx2, canony2;
	u64 xmin, ymin, xmax, ymax;

	canonicalize(ent->pos, &canonx1, &canony1);
	canonicalize(newpos, &canonx2, &canony2);

	xmin = MIN(canonx1, canonx2) - 1;
	ymin = MIN(canony1, canony2) - 1;
	xmax = MAX(canonx1, canonx2) + 1;
	ymax = MAX(canony1, canony2) + 1;

	u64 tiles[2][64] = {0};
	u64 dists[64] = {0};
	int n = 0;

	for(u64 i = ymin; i <= ymax; ++i) {
		for(u64 j = xmin; j <= xmax; ++j) {
			if(CHUNK[i*CHUNK_WIDTH + j] == 1) {
				tiles[1][n] = i;
				tiles[0][n] = j;
				dists[n] = abs(j - canonx1) + abs(i - canony1);
				++n;
			}
		}
	}

	// sort tiles by distance
	for(int i = 0; i < n - 1; ++i) {
		for(int j = 0; j < n - i - 1; ++j) {
			if(dists[j] > dists[j + 1]) {
				u64 tmp = dists[j];
				dists[j] = dists[j+1];
				dists[j+1] = tmp;
				tmp = tiles[0][j];
				tiles[0][j] = tiles[0][j+1];
				tiles[0][j+1] = tmp;
				tmp = tiles[1][j];
				tiles[1][j] = tiles[1][j+1];
				tiles[1][j+1] = tmp;
			}
		}
	}

	int pencount = 2;
	Vector2 a = ent->pos;
	//float t = Vector2Length(offset);
	float t = timestep;
	float invt = 1/t;
	Vector2 dir = { offset.x*invt, offset.y*invt };
	Vector2 invdir = { t/offset.x, t/offset.y };
	for(int i = 0; i < n; ++i) {
		float tx0 = (float)(tiles[0][i]*TILE_SIZE) - ent->widthH;
		float ty0 = (float)(tiles[1][i]*TILE_SIZE) - ent->heightH;
		float twidth = DOBL(ent->widthH) + TILE_SIZE;
		float theight = DOBL(ent->heightH) + TILE_SIZE;
		float tx1 = tx0 + twidth;
		float ty1 = ty0 + theight;

		Vector2 min = { tx0, ty0 }, max = { tx1, ty1 };

		float toi = 0;
		bool collision = rayAABBIntersect(ent->pos, invdir, min, max, &toi);

		if(collision && toi <= t) {
			Vector2 point = Vector2Add(ent->pos, Vector2Scale(dir, toi));
			Vector2 c = { tx0 + ent->widthH + HALF_TILE_SIZE, ty0 + ent->heightH + HALF_TILE_SIZE };
			Vector2 normal = Vector2Subtract(point, c);
			if(fabsf(normal.x) < fabsf(normal.y)) {
				newpos.y = point.y;
				float dy = newpos.y - a.y;
				dir.y = dy*invt;
				invdir.y = t/dy;
			} else {
				newpos.x = point.x;
				float dx = newpos.x - a.x;
				dir.x = dx*invt;
				invdir.x = t/dx;
			}
			CHUNK[tiles[1][i]*CHUNK_WIDTH + tiles[0][i]] = pencount++;
		}
	}

	float levelxmax = TILE_SIZE*CHUNK_WIDTH + TILE_SIZE;
	float levelymax = TILE_SIZE*CHUNK_HEIGHT + TILE_SIZE;

	if(newpos.x - ent->widthH < 0) {
		*(u32*)&ent->vel.x = NEGATIVE_ZERO;
		newpos.x = ent->widthH;
	} else if(newpos.x + ent->widthH > levelxmax) {
		ent->vel.x = 0.0;
		newpos.x = levelxmax - ent->widthH;
	}

	if(newpos.y - ent->heightH < 0) {
		ent->vel.y = 0;
		newpos.y = ent->heightH;
	} else if(newpos.y + ent->heightH > levelymax) {
		ent->vel.y = 0;
		newpos.y = levelymax - ent->heightH;
	}

	ent->pos = newpos;
}

INLINE void entityDraw(Entity *ent) {
	Rectangle entrec = {
		ent->pos.x - ent->widthH, ent->pos.y - ent->heightH, DOBL(ent->widthH), DOBL(ent->heightH),
	};

	DrawRectangleRec(entrec, RED);
}

int main(void) {
	char buf[128] = {0};
	SetWindowState(FLAG_WINDOW_RESIZABLE);
	InitWindow(1000, 800, "Tile Engine");
	SetTargetFPS(60);

	Entity bob = {
		.pos = { 200, 200 },
		.invmass = 1.0/27.0,
		.widthH = HALF(TILE_SIZE),
		.heightH = HALF(TILE_SIZE),
	};

	for(; !WindowShouldClose(); EndDrawing()) {
		BeginDrawing();

		float timestep = GetFrameTime();

		ClearBackground(BLACK);

		/* get input */
		entityGetInput(&bob);

		/* simulate */
		entitySimulate(&bob, timestep);

		/* draw */
		entityDraw(&bob);

		DrawLineEx(bob.pos, Vector2Add(bob.pos, bob.accel), 2.0, YELLOW);

		// draw tiles
		Rectangle tile = { 0, 0, TILE_SIZE, TILE_SIZE };
		for(int i = 0; i < 32; ++i) {
			for(int j = 0; j < 32; ++j) {
				int idx = i*CHUNK_WIDTH + j;
				if(CHUNK[idx] > 0) {
					tile.x = j*TILE_SIZE;
					tile.y = i*TILE_SIZE;
					if(CHUNK[idx] == 2) {
						DrawRectangleLinesEx(tile, 2.0, BLUE);
						CHUNK[idx] = 1;
					} else if(CHUNK[idx] > 2) {
						DrawRectangleLinesEx(tile, 2.0, GREEN);
						CHUNK[idx] = 1;
					} else {
						DrawRectangleLinesEx(tile, 2.0, RAYWHITE);
					}
				}
			}
		}

		stbsp_sprintf(buf, "FPS: %i\nFRAMETIME: %f\nforce.x: %f\nforce.y: %f\naccel.x: %f\naccel.y: %f\nvel.x: %f\nvel.y: %f\npos.x: %f\npos.y: %f\n",
				GetFPS(), timestep, bob.force.x, bob.force.y, bob.accel.x, bob.accel.y, bob.vel.x, bob.vel.y, bob.pos.x, bob.pos.y);
		DrawText(buf, 10, 10, 17, RAYWHITE);
	}

	CloseWindow();
	return 0;
}
