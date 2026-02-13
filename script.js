const canvas = document.getElementById("game");
const ctx = canvas.getContext("2d");
const scoreEl = document.getElementById("score");
const highScoreEl = document.getElementById("high-score");
const messageEl = document.getElementById("message");

const gridSize = 21;
const tileSize = canvas.width / gridSize;
const baseSpeed = 140;

let snake;
let direction;
let queuedDirection;
let food;
let score;
let highScore = Number(localStorage.getItem("snake-high-score")) || 0;
let gameOver;
let started;
let moveDelay;
let lastStep = 0;

highScoreEl.textContent = String(highScore);
resetGame();
requestAnimationFrame(gameLoop);

window.addEventListener("keydown", (event) => {
  const next = keyToDirection(event.key);
  if (next) {
    if (!started) {
      started = true;
      messageEl.textContent = "";
    }

    if (!isOpposite(next, direction)) {
      queuedDirection = next;
    }
  }

  if (event.code === "Space" && gameOver) {
    resetGame();
  }
});

function keyToDirection(key) {
  switch (key) {
    case "ArrowUp":
    case "w":
    case "W":
      return { x: 0, y: -1 };
    case "ArrowDown":
    case "s":
    case "S":
      return { x: 0, y: 1 };
    case "ArrowLeft":
    case "a":
    case "A":
      return { x: -1, y: 0 };
    case "ArrowRight":
    case "d":
    case "D":
      return { x: 1, y: 0 };
    default:
      return null;
  }
}

function isOpposite(a, b) {
  return a.x === -b.x && a.y === -b.y;
}

function resetGame() {
  snake = [
    { x: 10, y: 10 },
    { x: 9, y: 10 },
    { x: 8, y: 10 }
  ];
  direction = { x: 1, y: 0 };
  queuedDirection = direction;
  food = spawnFood();
  score = 0;
  gameOver = false;
  started = false;
  moveDelay = baseSpeed;
  lastStep = 0;
  scoreEl.textContent = "0";
  messageEl.textContent = "Press any arrow key to start.";
  draw();
}

function spawnFood() {
  while (true) {
    const candidate = {
      x: Math.floor(Math.random() * gridSize),
      y: Math.floor(Math.random() * gridSize)
    };

    if (!snake?.some((segment) => segment.x === candidate.x && segment.y === candidate.y)) {
      return candidate;
    }
  }
}

function gameLoop(timestamp) {
  requestAnimationFrame(gameLoop);

  if (!started || gameOver) {
    return;
  }

  if (timestamp - lastStep < moveDelay) {
    return;
  }

  lastStep = timestamp;
  step();
  draw();
}

function step() {
  if (!isOpposite(queuedDirection, direction)) {
    direction = queuedDirection;
  }

  const head = snake[0];
  const next = {
    x: (head.x + direction.x + gridSize) % gridSize,
    y: (head.y + direction.y + gridSize) % gridSize
  };

  const hitSelf = snake.some((segment) => segment.x === next.x && segment.y === next.y);
  if (hitSelf) {
    gameOver = true;
    messageEl.textContent = "Game over! Press Space to restart.";
    return;
  }

  snake.unshift(next);

  if (next.x === food.x && next.y === food.y) {
    score += 1;
    scoreEl.textContent = String(score);
    moveDelay = Math.max(70, baseSpeed - score * 4);
    food = spawnFood();
    if (score > highScore) {
      highScore = score;
      localStorage.setItem("snake-high-score", String(highScore));
      highScoreEl.textContent = String(highScore);
    }
  } else {
    snake.pop();
  }
}

function draw() {
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  drawGrid();

  ctx.fillStyle = "#ef4444";
  drawTile(food.x, food.y);

  snake.forEach((segment, index) => {
    ctx.fillStyle = index === 0 ? "#22c55e" : "#16a34a";
    drawTile(segment.x, segment.y);
  });

  if (gameOver) {
    ctx.fillStyle = "#000a";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = "#f8fafc";
    ctx.font = "bold 36px sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Game Over", canvas.width / 2, canvas.height / 2);
  }
}

function drawGrid() {
  ctx.strokeStyle = "#1e293b";
  ctx.lineWidth = 1;

  for (let i = 0; i <= gridSize; i += 1) {
    const p = i * tileSize;
    ctx.beginPath();
    ctx.moveTo(p, 0);
    ctx.lineTo(p, canvas.height);
    ctx.stroke();

    ctx.beginPath();
    ctx.moveTo(0, p);
    ctx.lineTo(canvas.width, p);
    ctx.stroke();
  }
}

function drawTile(x, y) {
  const padding = 1.5;
  ctx.fillRect(
    x * tileSize + padding,
    y * tileSize + padding,
    tileSize - padding * 2,
    tileSize - padding * 2
  );
}
