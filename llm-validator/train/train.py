import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

# Step 1: Load Data
# Assume you have these arrays ready
X_gene = np.load('X_gene.npy')   # shape (n_samples, 5000)
X_drug = np.load('X_drug.npy')   # shape (n_samples, 1024)
y_ic50 = np.load('y_ic50.npy')   # shape (n_samples,)

# Combine features
X_full = np.concatenate([X_gene, X_drug], axis=1)  # (n_samples, 6024)

# Convert to PyTorch tensors
X_tensor = torch.tensor(X_full, dtype=torch.float32)
y_tensor = torch.tensor(y_ic50, dtype=torch.float32).unsqueeze(1)

# Create Dataset and DataLoader
dataset = TensorDataset(X_tensor, y_tensor)
train_loader = DataLoader(dataset, batch_size=32, shuffle=True)

# Step 2: Define Model
class DrugResponseModel(nn.Module):
    def __init__(self, input_dim):
        super(DrugResponseModel, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, 2048),
            nn.ReLU(),
            nn.Linear(2048, 512),
            nn.ReLU(),
            nn.Linear(512, 1)  # Output: IC50 (regression)
        )

    def forward(self, x):
        return self.model(x)

input_dim = X_full.shape[1]  # 5000 + 1024 = 6024
model = DrugResponseModel(input_dim)

# Step 3: Define Loss and Optimizer
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=1e-4)

# Step 4: Training Loop
n_epochs = 50

for epoch in range(n_epochs):
    model.train()
    running_loss = 0.0
    for X_batch, y_batch in train_loader:
        optimizer.zero_grad()
        outputs = model(X_batch)
        loss = criterion(outputs, y_batch)
        loss.backward()
        optimizer.step()
        running_loss += loss.item()
    
    epoch_loss = running_loss / len(train_loader)
    print(f"Epoch {epoch+1}/{n_epochs}, Loss: {epoch_loss:.4f}")

print("Training finished!")
