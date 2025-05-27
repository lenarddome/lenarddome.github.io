---
layout: about
title: home
permalink: /
description: computational cognitive scientist
cv_pdf: lenarddome_cv.pdf
# profile:  
#    align: right
#    image: banner.png
#    address: lenarddome[at]gmail[dot]com
news: true  # includes a list of news items
selected_papers: true # includes a list of papers marked as selected={true}
social: false  # includes social icons at the bottom of the page
---
<br>


<div align="center" style="display: flex; flex-wrap: wrap; justify-content: center; gap: 32px; margin-bottom: 24px;">
   <div style="display: flex; flex-direction: column; align-items: center;">
      <div style="font-weight:bold; margin-bottom:4px;">Game of Life</div>
      <canvas id="gameOfLife" width="300" height="300" style="border:1px solid #ccc;"></canvas>
   </div>
   <div style="display: flex; flex-direction: column; align-items: center;">
      <div style="font-weight:bold; margin-bottom:4px;">Boids Swarm</div>
      <canvas id="swarmCanvas" width="300" height="300" style="border:1px solid #ccc;"></canvas>
   </div>
   <div style="display: flex; flex-direction: column; align-items: center;">
      <div style="font-weight:bold; margin-bottom:4px;">Neural Network</div>
      <canvas id="nnCanvas" width="300" height="300" style="border:1px solid #ccc; margin-top:0;"></canvas>
   </div>
</div>

{% raw %}
<script>
const canvas = document.getElementById('gameOfLife');
const ctx = canvas.getContext('2d');
const size = 40; // 40x40 grid
const cellSize = canvas.width / size;
let grid = Array.from({length: size}, () => Array(size).fill(0));

// Initialize with random cells
for (let y = 0; y < size; y++) {
   for (let x = 0; x < size; x++) {
      grid[y][x] = Math.random() > 0.7 ? 1 : 0;
   }
}

function draw() {
   ctx.clearRect(0, 0, canvas.width, canvas.height);
   for (let y = 0; y < size; y++) {
      for (let x = 0; x < size; x++) {
         ctx.fillStyle = grid[y][x] ? '#222' : '#fff';
         ctx.fillRect(x * cellSize, y * cellSize, cellSize, cellSize);
      }
   }
}

function nextGen() {
   const next = Array.from({length: size}, () => Array(size).fill(0));
   for (let y = 0; y < size; y++) {
      for (let x = 0; x < size; x++) {
         let neighbors = 0;
         for (let dy = -1; dy <= 1; dy++) {
            for (let dx = -1; dx <= 1; dx++) {
               if (dx === 0 && dy === 0) continue;
               const ny = (y + dy + size) % size;
               const nx = (x + dx + size) % size;
               neighbors += grid[ny][nx];
            }
         }
         if (grid[y][x]) {
            next[y][x] = neighbors === 2 || neighbors === 3 ? 1 : 0;
         } else {
            next[y][x] = neighbors === 3 ? 1 : 0;
         }
      }
   }
   grid = next;
}

function loop() {
   draw();
   nextGen();
   requestAnimationFrame(loop);
}

loop();
</script>
{% endraw %}

<div align="center" style="margin-top: 20px;">
</div>

{% raw %}
<script>
document.addEventListener('DOMContentLoaded', function() {
   const swarmCanvas = document.getElementById('swarmCanvas');
   const swarmCtx = swarmCanvas.getContext('2d');
   const swarmSize = 50;
   const boids = [];

   function randomVec() {
      const angle = Math.random() * 2 * Math.PI;
      return {x: Math.cos(angle), y: Math.sin(angle)};
   }

   for (let i = 0; i < swarmSize; i++) {
      boids.push({
         x: Math.random() * swarmCanvas.width,
         y: Math.random() * swarmCanvas.height,
         vx: (Math.random() - 0.5) * 2,
         vy: (Math.random() - 0.5) * 2
      });
   }

   function distance(a, b) {
      return Math.hypot(a.x - b.x, a.y - b.y);
   }

   function updateBoids() {
      const alignDist = 40, cohDist = 40, sepDist = 20;
      const alignWeight = 1, cohWeight = 0.5, sepWeight = 1.5, maxSpeed = 2;

      for (let i = 0; i < boids.length; i++) {
         let align = {x: 0, y: 0}, coh = {x: 0, y: 0}, sep = {x: 0, y: 0};
         let totalAlign = 0, totalCoh = 0, totalSep = 0;
         const b = boids[i];

         for (let j = 0; j < boids.length; j++) {
            if (i === j) continue;
            const other = boids[j];
            const d = distance(b, other);

            if (d < alignDist) {
               align.x += other.vx;
               align.y += other.vy;
               totalAlign++;
            }
            if (d < cohDist) {
               coh.x += other.x;
               coh.y += other.y;
               totalCoh++;
            }
            if (d < sepDist) {
               sep.x += b.x - other.x;
               sep.y += b.y - other.y;
               totalSep++;
            }
         }

         if (totalAlign) {
            align.x /= totalAlign; align.y /= totalAlign;
            const mag = Math.hypot(align.x, align.y) || 1;
            align.x = (align.x / mag) * alignWeight;
            align.y = (align.y / mag) * alignWeight;
         }
         if (totalCoh) {
            coh.x = (coh.x / totalCoh - b.x) * cohWeight / 100;
            coh.y = (coh.y / totalCoh - b.y) * cohWeight / 100;
         }
         if (totalSep) {
            sep.x = (sep.x / totalSep) * sepWeight;
            sep.y = (sep.y / totalSep) * sepWeight;
         }

         b.vx += align.x + coh.x + sep.x;
         b.vy += align.y + coh.y + sep.y;

         // Limit speed
         const speed = Math.hypot(b.vx, b.vy);
         if (speed > maxSpeed) {
            b.vx = (b.vx / speed) * maxSpeed;
            b.vy = (b.vy / speed) * maxSpeed;
         }
      }

      // Move boids and wrap around edges
      for (const b of boids) {
         b.x += b.vx;
         b.y += b.vy;
         if (b.x < 0) b.x += swarmCanvas.width;
         if (b.x > swarmCanvas.width) b.x -= swarmCanvas.width;
         if (b.y < 0) b.y += swarmCanvas.height;
         if (b.y > swarmCanvas.height) b.y -= swarmCanvas.height;
      }
   }

   function drawBoids() {
      swarmCtx.clearRect(0, 0, swarmCanvas.width, swarmCanvas.height);
      for (const b of boids) {
         swarmCtx.save();
         swarmCtx.translate(b.x, b.y);
         swarmCtx.rotate(Math.atan2(b.vy, b.vx));
         swarmCtx.beginPath();
         swarmCtx.moveTo(8, 0);
         swarmCtx.lineTo(-6, 4);
         swarmCtx.lineTo(-6, -4);
         swarmCtx.closePath();
         swarmCtx.fillStyle = "#0074D9";
         swarmCtx.fill();
         swarmCtx.restore();
      }
   }

   function swarmLoop() {
      updateBoids();
      drawBoids();
      requestAnimationFrame(swarmLoop);
   }

   swarmLoop();
});
</script>
{% endraw %}

{% raw %}
<script>
// Neural network with pulsing (spiking) neurons
const nnCanvas = document.getElementById('nnCanvas');
const nnCtx = nnCanvas.getContext('2d');
const neuronCount = 12;
const neurons = [];
const connections = [];
const radius = 110; // Reduced radius to fit 300x300
const center = {x: nnCanvas.width/2, y: nnCanvas.height/2};
const pulseDuration = 10; // frames

// Create neurons in a circle
for (let i = 0; i < neuronCount; i++) {
  const angle = (i / neuronCount) * 2 * Math.PI;
  neurons.push({
    x: center.x + radius * Math.cos(angle),
    y: center.y + radius * Math.sin(angle),
    spiking: false,
    pulseTimer: 0,
    id: i
  });
}

// Randomly connect neurons (each neuron connects to 2-3 others)
for (let i = 0; i < neuronCount; i++) {
  let targets = [];
  while (targets.length < 3) {
    let t = Math.floor(Math.random() * neuronCount);
    if (t !== i && !targets.includes(t)) targets.push(t);
  }
  for (let t of targets) {
    // Randomly assign excitatory (true) or inhibitory (false)
    const excitatory = Math.random() < 0.7; // 70% excitatory, 30% inhibitory
    connections.push({from: i, to: t, pulse: 0, excitatory});
  }
}

// Add inhibition state to neurons
for (const n of neurons) {
  n.inhibited = 0; // frames remaining inhibited
}

function drawNeuralNet() {
  nnCtx.clearRect(0, 0, nnCanvas.width, nnCanvas.height);
  // Draw connections
  for (const conn of connections) {
    const from = neurons[conn.from];
    const to = neurons[conn.to];
    if (conn.excitatory) {
      nnCtx.strokeStyle = conn.pulse > 0 ? '#00C853' : '#2196F3'; // green/blue for excitation
    } else {
      nnCtx.strokeStyle = conn.pulse > 0 ? '#FF1744' : '#B71C1C'; // red for inhibition
    }
    nnCtx.lineWidth = conn.pulse > 0 ? 3 : 1;
    nnCtx.beginPath();
    nnCtx.moveTo(from.x, from.y);
    nnCtx.lineTo(to.x, to.y);
    nnCtx.stroke();
  }
  // Draw neurons
  for (const n of neurons) {
    nnCtx.beginPath();
    nnCtx.arc(n.x, n.y, 12, 0, 2 * Math.PI);
    if (n.inhibited > 0) {
      nnCtx.fillStyle = '#888'; // gray for inhibited
    } else {
      nnCtx.fillStyle = n.spiking ? '#FFDC00' : '#0074D9';
    }
    nnCtx.shadowColor = n.spiking ? '#FFDC00' : 'transparent';
    nnCtx.shadowBlur = n.spiking ? 20 : 0;
    nnCtx.fill();
    nnCtx.shadowBlur = 0;
    nnCtx.strokeStyle = '#222';
    nnCtx.lineWidth = 2;
    nnCtx.stroke();
  }
}

function updateNeuralNet() {
  // Randomly spike neurons (only if not inhibited)
  for (const n of neurons) {
    if (n.inhibited > 0) {
      n.inhibited--;
      n.spiking = false;
      continue;
    }
    if (!n.spiking && Math.random() < 0.02) {
      n.spiking = true;
      n.pulseTimer = pulseDuration;
      // Send pulse to connections
      for (const conn of connections) {
        if (conn.from === n.id) conn.pulse = pulseDuration;
      }
    }
    if (n.spiking) {
      n.pulseTimer--;
      if (n.pulseTimer <= 0) n.spiking = false;
    }
  }
  // Propagate pulses
  for (const conn of connections) {
    if (conn.pulse > 0) {
      conn.pulse--;
      if (conn.pulse === Math.floor(pulseDuration/2)) {
        const target = neurons[conn.to];
        if (conn.excitatory) {
          // Excitatory: spike if not inhibited
          if (!target.spiking && target.inhibited === 0) {
            target.spiking = true;
            target.pulseTimer = pulseDuration;
          }
        } else {
          // Inhibitory: suppress target
          target.spiking = false;
          target.inhibited = pulseDuration * 2; // Inhibit for longer
        }
      }
    }
  }
}

function nnLoop() {
  updateNeuralNet();
  drawNeuralNet();
  requestAnimationFrame(nnLoop);
}

nnLoop();
</script>
{% endraw %}
